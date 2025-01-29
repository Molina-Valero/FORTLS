#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <limits>  // For NaN
#include <map>
#include <string>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace Rcpp;
using namespace std;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// euclidean distance function
inline double getEuclideanDistance(double x1, double y1, double z1, 
                                   double x2, double y2, double z2) {
  double dx = x2 - x1;
  double dy = y2 - y1;
  double dz = z2 - z1;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}


// geometric_features_point function

// [[Rcpp::export]]
std::map<std::string, double> geometric_features_point(const Eigen::Map<Eigen::VectorXd>& x, 
                                                       const Eigen::Map<Eigen::VectorXd>& y, 
                                                       const Eigen::Map<Eigen::VectorXd>& z, 
                                                       int pto, 
                                                       double dist,
                                                       bool First_eigenvalue = true,
                                                       bool Second_eigenvalue = true,
                                                       bool Third_eigenvalue = true,
                                                       bool Sum_of_eigenvalues = true,
                                                       bool PCA_1 = true,
                                                       bool PCA_2 = true,
                                                       bool Anisotropy = true,
                                                       bool Planarity = true,
                                                       bool Linearity = true,
                                                       bool Surface_variation = true,
                                                       bool Normal_change_rate = true,
                                                       bool Verticality = true,
                                                       bool Number_of_points = true,
                                                       bool omnivariance = true,
                                                       bool eigenentropy = true,
                                                       bool surface_density = true,
                                                       bool volume_density = true,
                                                       int solver_thresh = 50000) {      // 'solver_thresh': number of chunk-points (threshold) where the solver will switch from 'compute' to 'computeDirect' (approximate)
  const int SIZE = x.size();
  
  // First pass: count points within distance
  int count = 0;
  std::vector<int> indices_in_range;
  indices_in_range.reserve(SIZE);  // Pre-allocate max possible size to avoid resizing
  
  // Collect indices of points within the distance from the point of interest
  for (int i = 0; i < SIZE; ++i) {
    double distance = getEuclideanDistance(x[pto], y[pto], z[pto], x[i], y[i], z[i]);
    if (distance <= dist) {
      indices_in_range.push_back(i);
      ++count;
    }
  }
  
  // Assign features to a map
  std::map<std::string, double> features;
  
  if (count == 0) {
    
    if (First_eigenvalue) features["First_eigenvalue"] = std::numeric_limits<double>::quiet_NaN();
    if (Second_eigenvalue) features["Second_eigenvalue"] = std::numeric_limits<double>::quiet_NaN();
    if (Third_eigenvalue) features["Third_eigenvalue"] = std::numeric_limits<double>::quiet_NaN();
    if (Sum_of_eigenvalues) features["Sum_of_eigenvalues"] = std::numeric_limits<double>::quiet_NaN();
    if (PCA_1) features["PCA1"] = std::numeric_limits<double>::quiet_NaN();
    if (PCA_2) features["PCA2"] = std::numeric_limits<double>::quiet_NaN();
    if (Anisotropy) features["Anisotropy"] = std::numeric_limits<double>::quiet_NaN();
    if (Planarity) features["Planarity"] = std::numeric_limits<double>::quiet_NaN();
    if (Linearity) features["Linearity"] = std::numeric_limits<double>::quiet_NaN();
    if (Surface_variation) features["Surface_variation"] = std::numeric_limits<double>::quiet_NaN();
    if (Normal_change_rate) features["Normal_change_rate"] = std::numeric_limits<double>::quiet_NaN();
    if (Verticality) features["Verticality"] = std::numeric_limits<double>::quiet_NaN();
    if (Number_of_points) features["Number_of_points"] = std::numeric_limits<double>::quiet_NaN();
    if (omnivariance) features["omnivariance"] = std::numeric_limits<double>::quiet_NaN();
    if (eigenentropy) features["eigenentropy"] = std::numeric_limits<double>::quiet_NaN();
    if (surface_density) features["surface_density"] = std::numeric_limits<double>::quiet_NaN();
    if (volume_density) features["volume_density"] = std::numeric_limits<double>::quiet_NaN();
    
    // Return an empty map if no points are within the distance
    return features;
  }
  
  // Allocate matrix based on the number of points within range
  Eigen::MatrixXd m(count, 3);
  
  // Fill the matrix with the coordinates of points within the range
  for (int idx = 0; idx < count; ++idx) {
    int i = indices_in_range[idx];
    m(idx, 0) = x[i];
    m(idx, 1) = y[i];
    m(idx, 2) = z[i];
  }
  
  // Calculate the covariance matrix
  Eigen::MatrixXd centered = m.rowwise() - m.colwise().mean();
  Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(m.rows() - 1);
  
  // Compute eigenvalues (no eigenvectors yet)
  Eigen::SelfAdjointEigenSolver<MatrixXd> es;
  // es.compute(cov, true);                            // Compute both eigenvalues and eigenvectors
  
  // Compute both eigenvalues and eigenvectors
  if (SIZE > solver_thresh) {
    es.computeDirect(cov, true);
  } else {
    es.compute(cov, true);
  }
  
  Eigen::VectorXd eigenvalues = es.eigenvalues().real();
  
  // Compute the verticality only if required
  double ver;
  if (Verticality) {

    Eigen::MatrixXd eigenvectors = es.eigenvectors().real();
    int minEigenIndex = 0;
    
    for (int i = 1; i < 3; ++i) {
      if (eigenvalues(i) < eigenvalues(minEigenIndex)) {
        minEigenIndex = i;
      }
    }
    
    Eigen::Vector3d normalVector = eigenvectors.col(minEigenIndex);
    normalVector.normalize();
    
    ver = 1.0 - std::abs(normalVector[2]);
  } else {
    ver = std::numeric_limits<double>::quiet_NaN();
  }
  
  // Precompute values for reuse
  double lambda1 = (First_eigenvalue) ? eigenvalues(2) : std::numeric_limits<double>::quiet_NaN();
  double lambda2 = (Second_eigenvalue) ? eigenvalues(1) : std::numeric_limits<double>::quiet_NaN();
  double lambda3 = (Third_eigenvalue) ? eigenvalues(0) : std::numeric_limits<double>::quiet_NaN();
  double lambda123 = (Sum_of_eigenvalues) ? eigenvalues.sum() : std::numeric_limits<double>::quiet_NaN();
  
  features["First_eigenvalue"] = lambda1;
  features["Second_eigenvalue"] = lambda2;
  features["Third_eigenvalue"] = lambda3;
  features["Sum_of_eigenvalues"] = lambda123;
  features["PCA1"] = (PCA_1) ? lambda1 / lambda123 : std::numeric_limits<double>::quiet_NaN();
  features["PCA2"] = (PCA_2) ? lambda2 / lambda123 : std::numeric_limits<double>::quiet_NaN();
  features["Anisotropy"] = (Anisotropy) ? (lambda1 - lambda3) / lambda1 : std::numeric_limits<double>::quiet_NaN();
  features["Planarity"] = (Planarity) ? (lambda2 - lambda3) / lambda1 : std::numeric_limits<double>::quiet_NaN();
  features["Linearity"] = (Linearity) ? (lambda1 - lambda2) / lambda1 : std::numeric_limits<double>::quiet_NaN();
  features["Surface_variation"] = (Surface_variation) ? lambda3 / lambda123 : std::numeric_limits<double>::quiet_NaN();
  features["Normal_change_rate"] = (Normal_change_rate) ? lambda3 / lambda1 : std::numeric_limits<double>::quiet_NaN();
  features["Verticality"] = ver;
  features["Number_of_points"] = (Number_of_points) ? static_cast<double>(count) : std::numeric_limits<double>::quiet_NaN();
  features["omnivariance"] = (omnivariance) ? std::pow(lambda1 * lambda2 * lambda3, 1.0 / 3.0)  : std::numeric_limits<double>::quiet_NaN();
  features["eigenentropy"] = (eigenentropy) ? -(lambda1 * std::log(lambda1) + lambda2 * std::log(lambda2) + lambda3 * std::log(lambda3)) : std::numeric_limits<double>::quiet_NaN();
  features["surface_density"] = (surface_density) ? static_cast<double>(count) / (4.0 * M_PI * std::pow(dist, 2.0))  : std::numeric_limits<double>::quiet_NaN();
  features["volume_density"] = (volume_density) ? static_cast<double>(count) / ((4.0/3.0) * M_PI * std::pow(dist, 3.0))  : std::numeric_limits<double>::quiet_NaN();

  return features;
}


// geometric_features function

// [[Rcpp::export]]
DataFrame geometric_features(const Eigen::MatrixXd& m,
                             double dist,
                             bool First_eigenvalue = true,
                             bool Second_eigenvalue = true,
                             bool Third_eigenvalue = true,
                             bool Sum_of_eigenvalues = true,
                             bool PCA_1 = true,
                             bool PCA_2 = true,
                             bool Anisotropy = true,
                             bool Planarity = true,
                             bool Linearity = true,
                             bool Surface_variation = true,
                             bool Normal_change_rate = true,
                             bool Verticality = true,
                             bool Number_of_points = true,
                             bool omnivariance = true,
                             bool eigenentropy = true,
                             bool surface_density = true,
                             bool volume_density = true,
                             int threads = 1,
                             int solver_thresh = 50000) {      // 'solver_thresh': number of chunk-points (threshold) where the solver will switch from 'compute' to 'computeDirect' (approximate)
  #ifdef _OPENMP
  omp_set_num_threads(threads);
  #endif

  Eigen::VectorXd point = m.col(0);
  Eigen::VectorXd x = m.col(1);
  Eigen::VectorXd y = m.col(2);
  Eigen::VectorXd z = m.col(3);
  
  // Create Eigen::Map for each vector
  const Eigen::Map<Eigen::VectorXd> x_map(x.data(), x.size());
  const Eigen::Map<Eigen::VectorXd> y_map(y.data(), y.size());
  const Eigen::Map<Eigen::VectorXd> z_map(z.data(), z.size());

  double firstX = x.minCoeff();
  double lastX = x.maxCoeff();
  double firstY = y.minCoeff();
  double lastY = y.maxCoeff();
  double firstZ = z.minCoeff();
  double lastZ = z.maxCoeff();

  // Get the size of the input data
  unsigned int n = x.size();

  // Initialize vectors with NaN values using std::numeric_limits<double>::quiet_NaN()
  std::vector<double> lambda1(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> lambda2(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> lambda3(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> lambda123(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> PCA1_vec(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> PCA2_vec(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> ani(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> pla(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> lin(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> sur(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> sph(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> ver(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> pts(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> kk(n, std::numeric_limits<double>::quiet_NaN());    // although "kk" will be an integer vector the only way to assign a NaN is to convert to a double vector
  std::vector<double> omn(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> eigent(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> surfd(n, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> vold(n, std::numeric_limits<double>::quiet_NaN());

  unsigned int q;

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static) shared(x_map, y_map, z_map, dist, n, First_eigenvalue, Second_eigenvalue, Third_eigenvalue, Sum_of_eigenvalues, PCA_1, PCA_2, Anisotropy, Planarity, Linearity, Surface_variation, Normal_change_rate, Verticality, Number_of_points, omnivariance, eigenentropy, surface_density, volume_density, firstZ, lastZ, firstY, lastY, firstX, lastX, lambda1, lambda2, lambda3, lambda123, PCA1_vec, PCA2_vec, ani, pla, lin, sur, sph, ver, pts, kk, point, eigent, omn, surfd, vold, solver_thresh) default(none) private(q)
  #endif
  for (q = 0; q < n; ++q) {
    // Boundary check to avoid edge effects
    if (!(x_map[q] > lastX - dist || x_map[q] < firstX + dist || y_map[q] > lastY - dist || y_map[q] < firstY + dist || z_map[q] > lastZ - dist || z_map[q] < firstZ + dist)) {

     // Compute geometric features for the current point q
     std::map<std::string, double> features = geometric_features_point(x_map, y_map, z_map, q, dist,
                                                                       First_eigenvalue, Second_eigenvalue, Third_eigenvalue, Sum_of_eigenvalues,
                                                                       PCA_1, PCA_2, Anisotropy, Planarity, Linearity,
                                                                       Surface_variation, Normal_change_rate, Verticality, Number_of_points,
                                                                       omnivariance, eigenentropy, surface_density, volume_density, solver_thresh);

      // Using index assignment with OpenMP atomic to avoid race conditions
      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      lambda1[q] = features["First_eigenvalue"];

      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      lambda2[q] = features["Second_eigenvalue"];

      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      lambda3[q] = features["Third_eigenvalue"];

      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      lambda123[q] = features["Sum_of_eigenvalues"];

      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      PCA1_vec[q] = features["PCA1"];

      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      PCA2_vec[q] = features["PCA2"];

      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      ani[q] = features["Anisotropy"];

      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      pla[q] = features["Planarity"];

      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      lin[q] = features["Linearity"];

      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      sur[q] = features["Surface_variation"];

      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      sph[q] = features["Normal_change_rate"];

      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      ver[q] = features["Verticality"];

      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      pts[q] = features["Number_of_points"];

      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      kk[q] = static_cast<double>(point[q]);
      
      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      omn[q] = features["omnivariance"];
      
      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      eigent[q] = features["eigenentropy"];
      
      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      surfd[q] = features["surface_density"];
      
      #ifdef _OPENMP
      #pragma omp atomic write
      #endif
      vold[q] = features["volume_density"];
    }
  }

  // Create and return DataFrame with all results
  return DataFrame::create(Named("point") = kk,
                           Named("first_eigenvalue") = lambda1,
                           Named("second_eigenvalue") = lambda2,
                           Named("third_eigenvalue") = lambda3,
                           Named("sum_eigenvalues") = lambda123,
                           Named("PCA1") = PCA1_vec,
                           Named("PCA2") = PCA2_vec,
                           Named("anisotropy") = ani,
                           Named("planarity") = pla,
                           Named("linearity") = lin,
                           Named("surface_variation") = sur,
                           Named("sphericity") = sph,
                           Named("verticality") = ver,
                           Named("number_neighbors") = pts,
                           Named("omnivariance") = omn,
                           Named("eigenentropy") = eigent,
                           Named("surface_density") = surfd,
                           Named("volume_density") = vold);
}

