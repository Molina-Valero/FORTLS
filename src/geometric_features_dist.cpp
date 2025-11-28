// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <limits>
#include <array>
#include <algorithm>
#include <unordered_map>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// ----- SPATIAL INDEXING STRUCTURES -----

struct GridCell {
  int x, y, z;
  
  bool operator==(const GridCell& other) const {
    return x == other.x && y == other.y && z == other.z;
  }
};

struct GridCellHash {
  std::size_t operator()(const GridCell& cell) const {
    return static_cast<std::size_t>(cell.x) * 73856093 ^
           static_cast<std::size_t>(cell.y) * 19349663 ^
           static_cast<std::size_t>(cell.z) * 83492791;
  }
};

class SpatialGrid {
private:
  double cell_size_;
  double inv_cell_size_;  // Precompute inverse for faster division
  std::unordered_map<GridCell, std::vector<int>, GridCellHash> grid_;
  
public:
  explicit SpatialGrid(double cell_size) 
    : cell_size_(cell_size), inv_cell_size_(1.0 / cell_size) {
    grid_.reserve(1024);  // Pre-reserve to reduce rehashing
  }
  
  inline GridCell GetCell(double x, double y, double z) const {
    return GridCell{
      static_cast<int>(std::floor(x * inv_cell_size_)),
      static_cast<int>(std::floor(y * inv_cell_size_)),
      static_cast<int>(std::floor(z * inv_cell_size_))
    };
  }
  
  void Insert(int point_idx, double x, double y, double z) {
    GridCell cell = GetCell(x, y, z);
    grid_[cell].push_back(point_idx);
  }
  
  void Reserve() {
    for (auto& pair : grid_) {
      pair.second.shrink_to_fit();
    }
  }
  
  void QueryNeighbors(double x, double y, double z, double radius,
                      std::vector<int>& result) const {
    result.clear();
    
    GridCell center = GetCell(x, y, z);
    int cell_range = static_cast<int>(std::ceil(radius * inv_cell_size_)) + 1;
    
    // Estimate and reserve space
    result.reserve(cell_range * cell_range * cell_range * 8);
    
    for (int dx = -cell_range; dx <= cell_range; ++dx) {
      for (int dy = -cell_range; dy <= cell_range; ++dy) {
        for (int dz = -cell_range; dz <= cell_range; ++dz) {
          GridCell query_cell{center.x + dx, center.y + dy, center.z + dz};
          auto it = grid_.find(query_cell);
          if (it != grid_.end()) {
            result.insert(result.end(), it->second.begin(), it->second.end());
          }
        }
      }
    }
  }
};

// ----- PROGRESS BAR UTILITY -----
class ProgressBar {
private:
  int total_;
  int bar_width_;
  int last_percent_;
  std::string cached_bar_;
  
public:
  ProgressBar(int total, int bar_width = 50) 
    : total_(total), bar_width_(bar_width), last_percent_(-1) {
    cached_bar_.reserve(bar_width + 50);
  }
  
  void update(int current) {
    int percent = static_cast<int>((current * 100.0) / total_);
    
    if (percent == last_percent_) return;
    last_percent_ = percent;
    
    int filled = static_cast<int>((current * bar_width_) / total_);
    
    cached_bar_.clear();
    cached_bar_ += "\r[";
    for (int i = 0; i < bar_width_; ++i) {
      if (i < filled) cached_bar_ += '=';
      else if (i == filled) cached_bar_ += '>';
      else cached_bar_ += ' ';
    }
    cached_bar_ += "] " + std::to_string(percent) + "% (" + 
                   std::to_string(current) + "/" + std::to_string(total_) + ")";
    
    Rcpp::Rcout << cached_bar_ << std::flush;
  }
  
  void finish() {
    Rcpp::Rcout << std::endl;
  }
};

// ----- Type aliases -----
template <typename T>
using Vec3 = Eigen::Matrix<T, 3, 1>;

template <typename T>
using Mat3 = Eigen::Matrix<T, 3, 3>;

template <typename real_t>
struct PCAResult {
  Vec3<real_t> val;
  Vec3<real_t> v0, v1, v2;
};

// Compact struct for hot path - better cache utilization
struct GeometricFeatures {
  double point_id;
  double lambda1, lambda2, lambda3, eigenvalue_sum;
  double normal_x, normal_y, normal_z;
  double pca1, pca2;
  double anisotropy, eigenentropy, linearity;
  double omnivariance, planarity, sphericity;
  double surface_variation, verticality;
  double n_points, surface_density, volume_density;
  int count;
  
  GeometricFeatures() : 
    point_id(NA_REAL), lambda1(NA_REAL), lambda2(NA_REAL), lambda3(NA_REAL),
    eigenvalue_sum(NA_REAL), normal_x(NA_REAL), normal_y(NA_REAL), normal_z(NA_REAL),
    pca1(NA_REAL), pca2(NA_REAL), anisotropy(NA_REAL), eigenentropy(NA_REAL),
    linearity(NA_REAL), omnivariance(NA_REAL), planarity(NA_REAL), sphericity(NA_REAL),
    surface_variation(NA_REAL), verticality(NA_REAL), 
    n_points(NA_REAL), surface_density(NA_REAL), volume_density(NA_REAL), count(0) {}
};

// ----- AUXILIARY FUNCTIONS -----

inline double SquaredDistance(double x1, double y1, double z1,
                              double x2, double y2, double z2) {
  const double dx = x2 - x1;
  const double dy = y2 - y1;
  const double dz = z2 - z1;
  return dx * dx + dy * dy + dz * dz;
}

// ----- Optimized Eigenvalue analysis -----
template <typename real_t>
static inline PCAResult<real_t>
eigenvalue_analysis_core(const std::vector<Vec3<real_t>>& points, 
                        const Vec3<real_t>& centroid)
{
  const size_t N = points.size();
  
  // Incremental covariance computation (more cache-friendly)
  Mat3<real_t> cov = Mat3<real_t>::Zero();
  
  for (size_t i = 0; i < N; ++i) {
    const Vec3<real_t> diff = points[i] - centroid;
    // Manual outer product is faster than diff * diff.transpose()
    cov(0, 0) += diff(0) * diff(0);
    cov(1, 1) += diff(1) * diff(1);
    cov(2, 2) += diff(2) * diff(2);
    cov(0, 1) += diff(0) * diff(1);
    cov(0, 2) += diff(0) * diff(2);
    cov(1, 2) += diff(1) * diff(2);
  }
  
  // Symmetrize
  cov(1, 0) = cov(0, 1);
  cov(2, 0) = cov(0, 2);
  cov(2, 1) = cov(1, 2);
  cov /= static_cast<real_t>(N - 1);
  
  Eigen::SelfAdjointEigenSolver<Mat3<real_t>> es(cov);
  
  const auto ev = es.eigenvalues();
  const auto evecs = es.eigenvectors();
  
  // Simplified sorting network for 3 elements
  std::array<int, 3> idx{2, 1, 0};
  if (ev(1) > ev(2)) std::swap(idx[0], idx[1]);
  if (ev(0) > ev(1)) std::swap(idx[1], idx[2]);
  if (ev(1) > ev(2)) std::swap(idx[0], idx[1]);
  
  PCAResult<real_t> out;
  out.val(0) = std::max<real_t>(ev(idx[0]), real_t(0));
  out.val(1) = std::max<real_t>(ev(idx[1]), real_t(0));
  out.val(2) = std::max<real_t>(ev(idx[2]), real_t(0));
  
  out.v0 = evecs.col(idx[0]);
  out.v1 = evecs.col(idx[1]);
  out.v2 = evecs.col(idx[2]);
  
  if (out.v2(2) < real_t(0)) out.v2 = -out.v2;
  if (out.v0.cross(out.v1).dot(out.v2) < real_t(0)) out.v1 = -out.v1;
  
  return out;
}

// ----- Optimized Core feature computation -----
GeometricFeatures geometric_features_core(
    const double* x_data,
    const double* y_data,
    const double* z_data,
    const SpatialGrid& spatial_grid,
    double x_pto, double y_pto, double z_pto,
    double point_id, double dist,
    double dist_sq,
    bool compute_density,
    std::vector<int>& candidates_buffer,
    std::vector<Vec3<double>>& points_buffer)
{
  GeometricFeatures features;
  features.point_id = point_id;
  
  // Reuse buffers to avoid allocations
  spatial_grid.QueryNeighbors(x_pto, y_pto, z_pto, dist, candidates_buffer);
  
  // Filter and build point cloud in one pass
  points_buffer.clear();
  Vec3<double> centroid(0, 0, 0);
  
  for (int idx : candidates_buffer) {
    const double d_sq = SquaredDistance(x_pto, y_pto, z_pto, 
                                       x_data[idx], y_data[idx], z_data[idx]);
    if (d_sq <= dist_sq) {
      Vec3<double> p(x_data[idx], y_data[idx], z_data[idx]);
      points_buffer.push_back(p);
      centroid += p;
    }
  }
  
  const int count = static_cast<int>(points_buffer.size());
  features.count = count;
  features.n_points = static_cast<double>(count);
  
  if (count < 3) return features;
  
  centroid /= static_cast<double>(count);
  
  // PCA computation
  const auto pca = eigenvalue_analysis_core(points_buffer, centroid);
  
  const double lambda1 = pca.val(0);
  const double lambda2 = pca.val(1);
  const double lambda3 = pca.val(2);
  const double lambda_sum = lambda1 + lambda2 + lambda3;
  constexpr double epsilon = 1e-10;
  
  features.lambda1 = lambda1;
  features.lambda2 = lambda2;
  features.lambda3 = lambda3;
  features.eigenvalue_sum = lambda_sum;
  
  features.normal_x = pca.v2(0);
  features.normal_y = pca.v2(1);
  features.normal_z = pca.v2(2);
  
  features.pca1 = pca.v0(2);
  features.pca2 = pca.v1(2);
  
  // Geometric features (compute all at once to maximize instruction-level parallelism)
  if (lambda1 > epsilon) {
    const double inv_lambda1 = 1.0 / lambda1;
    features.linearity = (lambda1 - lambda2) * inv_lambda1;
    features.planarity = (lambda2 - lambda3) * inv_lambda1;
    features.sphericity = lambda3 * inv_lambda1;
    features.anisotropy = (lambda1 - lambda3) * inv_lambda1;
  }
  
  if (lambda_sum > epsilon) {
    features.surface_variation = lambda3 / lambda_sum;
  }
  
  if (lambda1 > epsilon && lambda2 > epsilon && lambda3 > epsilon) {
    const double inv_sum = 1.0 / lambda_sum;
    const double norm_l1 = lambda1 * inv_sum;
    const double norm_l2 = lambda2 * inv_sum;
    const double norm_l3 = lambda3 * inv_sum;
    features.eigenentropy = -(norm_l1 * std::log(norm_l1) + 
                              norm_l2 * std::log(norm_l2) + 
                              norm_l3 * std::log(norm_l3));
    features.omnivariance = std::cbrt(lambda1 * lambda2 * lambda3);
  }
  
  features.verticality = 1.0 - std::abs(pca.v2(2));
  
  if (compute_density && dist > epsilon) {
    constexpr double PI = 3.14159265358979323846;
    constexpr double FOUR_PI = 4.0 * PI;
    constexpr double FOUR_THIRDS_PI = (4.0 / 3.0) * PI;
    
    const double count_d = static_cast<double>(count);
    features.surface_density = count_d / (FOUR_PI * dist_sq);
    features.volume_density = count_d / (FOUR_THIRDS_PI * dist_sq * dist);
  }
  
  return features;
}

// [[Rcpp::export]]
Rcpp::DataFrame geometric_features_dist(
    Rcpp::NumericMatrix points,
    double dist,
    bool Anisotropy = false,
    bool Eigenentropy = false,
    bool Eigenvalue_sum = false,
    bool First_eigenvalue = false,
    bool Linearity = false,
    bool Normal_x = false,
    bool Normal_y = false,
    bool Normal_z = false,
    bool Number_of_points = false,
    bool Omnivariance = false,
    bool PCA_1 = false,
    bool PCA_2 = false,
    bool Planarity = false,
    bool Second_eigenvalue = false,
    bool Sphericity = false,
    bool Surface_density = false,
    bool Surface_variation = false,
    bool Third_eigenvalue = false,
    bool Verticality = false,
    bool Volume_density = false,
    int num_threads = 0,
    bool use_spatial_index = true,
    double grid_cell_size = -1.0) {
  
  if (points.ncol() != 4) {
    Rcpp::stop("Input 'points' must be a matrix with 4 columns: point, x, y, z");
  }
  
  const int n_points = points.nrow();
  if (n_points == 0) return Rcpp::DataFrame();
  
#ifdef _OPENMP
  if (num_threads > 0) {
    omp_set_num_threads(num_threads);
  }
#endif
  
  // Extract coordinates (more cache-friendly layout)
  std::vector<double> x_all(n_points), y_all(n_points), z_all(n_points);
  
  for (int i = 0; i < n_points; ++i) {
    x_all[i] = points(i, 1);
    y_all[i] = points(i, 2);
    z_all[i] = points(i, 3);
  }
  
  const double* x_data = x_all.data();
  const double* y_data = y_all.data();
  const double* z_data = z_all.data();
  
  const bool compute_density = Surface_density || Volume_density;
  const double dist_sq = dist * dist;
  
  // Build spatial index
  std::unique_ptr<SpatialGrid> spatial_grid;
  if (use_spatial_index) {
    double cell_size = grid_cell_size > 0 ? grid_cell_size : dist * 2.0;
    spatial_grid = std::make_unique<SpatialGrid>(cell_size);
    
    for (int i = 0; i < n_points; ++i) {
      spatial_grid->Insert(i, x_all[i], y_all[i], z_all[i]);
    }
    spatial_grid->Reserve();  // Optimize memory
  }
  
  // Pre-allocate result vectors
  std::vector<double> vec_point(n_points), vec_x(n_points), vec_y(n_points), vec_z(n_points);
  std::vector<double> vec_first_ev(n_points, NA_REAL), vec_second_ev(n_points, NA_REAL);
  std::vector<double> vec_third_ev(n_points, NA_REAL), vec_ev_sum(n_points, NA_REAL);
  std::vector<double> vec_normal_x(n_points, NA_REAL), vec_normal_y(n_points, NA_REAL);
  std::vector<double> vec_normal_z(n_points, NA_REAL);
  std::vector<double> vec_pca1(n_points, NA_REAL), vec_pca2(n_points, NA_REAL);
  std::vector<double> vec_anisotropy(n_points, NA_REAL), vec_eigenentropy(n_points, NA_REAL);
  std::vector<double> vec_linearity(n_points, NA_REAL), vec_omnivariance(n_points, NA_REAL);
  std::vector<double> vec_planarity(n_points, NA_REAL), vec_sphericity(n_points, NA_REAL);
  std::vector<double> vec_surface_var(n_points, NA_REAL);
  std::vector<double> vec_verticality(n_points, NA_REAL);
  std::vector<double> vec_n_points(n_points, NA_REAL);
  std::vector<double> vec_surf_density(n_points, NA_REAL), vec_vol_density(n_points, NA_REAL);
  
  // Process points in parallel
  Rcpp::Rcout << "Computing geometric features..." << std::endl;
  ProgressBar compute_progress(n_points);
  
#pragma omp parallel if(n_points > 100)
  {
    // Thread-local buffers to avoid allocations in hot loop
    std::vector<int> candidates_buffer;
    std::vector<Vec3<double>> points_buffer;
    candidates_buffer.reserve(500);
    points_buffer.reserve(500);
    
#pragma omp for schedule(dynamic, 64)
    for (int i = 0; i < n_points; ++i) {
      
      // Update progress bar (only from master thread)
#ifdef _OPENMP
      if (omp_get_thread_num() == 0 && i % 1000 == 0) {
        compute_progress.update(i);
      }
#else
      if (i % 1000 == 0) {
        compute_progress.update(i);
        Rcpp::checkUserInterrupt();
      }
#endif
      
      const double point_id = points(i, 0);
      const double x_pto = points(i, 1);
      const double y_pto = points(i, 2);
      const double z_pto = points(i, 3);
      
      // Compute features
      const auto features = geometric_features_core(
        x_data, y_data, z_data, *spatial_grid,
        x_pto, y_pto, z_pto, point_id, dist, dist_sq,
        compute_density, candidates_buffer, points_buffer);
      
      // Store results
      vec_point[i] = features.point_id;
      vec_x[i] = x_pto;
      vec_y[i] = y_pto;
      vec_z[i] = z_pto;
      
      if (First_eigenvalue) vec_first_ev[i] = features.lambda1;
      if (Second_eigenvalue) vec_second_ev[i] = features.lambda2;
      if (Third_eigenvalue) vec_third_ev[i] = features.lambda3;
      if (Eigenvalue_sum) vec_ev_sum[i] = features.eigenvalue_sum;
      if (Normal_x) vec_normal_x[i] = features.normal_x;
      if (Normal_y) vec_normal_y[i] = features.normal_y;
      if (Normal_z) vec_normal_z[i] = features.normal_z;
      if (PCA_1) vec_pca1[i] = features.pca1;
      if (PCA_2) vec_pca2[i] = features.pca2;
      if (Anisotropy) vec_anisotropy[i] = features.anisotropy;
      if (Eigenentropy) vec_eigenentropy[i] = features.eigenentropy;
      if (Linearity) vec_linearity[i] = features.linearity;
      if (Omnivariance) vec_omnivariance[i] = features.omnivariance;
      if (Planarity) vec_planarity[i] = features.planarity;
      if (Sphericity) vec_sphericity[i] = features.sphericity;
      if (Surface_variation) vec_surface_var[i] = features.surface_variation;
      if (Verticality) vec_verticality[i] = features.verticality;
      if (Number_of_points) vec_n_points[i] = features.n_points;
      if (Surface_density) vec_surf_density[i] = features.surface_density;
      if (Volume_density) vec_vol_density[i] = features.volume_density;
    }
  }
  
  compute_progress.update(n_points);
  compute_progress.finish();
  
  // Build DataFrame
  Rcpp::List result;
  result["point"] = vec_point;
  result["x"] = vec_x;
  result["y"] = vec_y;
  result["z"] = vec_z;
  
  if (Anisotropy) result["Anisotropy"] = vec_anisotropy;
  if (Eigenentropy) result["Eigenentropy"] = vec_eigenentropy;
  if (Eigenvalue_sum) result["Eigenvalue_sum"] = vec_ev_sum;
  if (First_eigenvalue) result["First_eigenvalue"] = vec_first_ev;
  if (Linearity) result["Linearity"] = vec_linearity;
  if (Normal_x) result["Normal_x"] = vec_normal_x;
  if (Normal_y) result["Normal_y"] = vec_normal_y;
  if (Normal_z) result["Normal_z"] = vec_normal_z;
  if (Number_of_points) result["Number_of_points"] = vec_n_points;
  if (Omnivariance) result["Omnivariance"] = vec_omnivariance;
  if (PCA_1) result["PCA_1"] = vec_pca1;
  if (PCA_2) result["PCA_2"] = vec_pca2;
  if (Planarity) result["Planarity"] = vec_planarity;
  if (Second_eigenvalue) result["Second_eigenvalue"] = vec_second_ev;
  if (Sphericity) result["Sphericity"] = vec_sphericity;
  if (Surface_density) result["Surface_density"] = vec_surf_density;
  if (Surface_variation) result["Surface_variation"] = vec_surface_var;
  if (Third_eigenvalue) result["Third_eigenvalue"] = vec_third_ev;
  if (Verticality) result["Verticality"] = vec_verticality;
  if (Volume_density) result["Volume_density"] = vec_vol_density;
  
  return Rcpp::DataFrame(result);
}
