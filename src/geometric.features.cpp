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

// Euclidean distance function
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
                                                       double x_pto,
                                                       double y_pto,
                                                       double z_pto,
                                                       double point_pto,
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
    double distance = getEuclideanDistance(x_pto, y_pto, z_pto, x[i], y[i], z[i]);
    // std::cout << distance << std::endl;

    if (distance <= dist) {
      indices_in_range.push_back(i);
      ++count;
    }
  }

  // Assign features to a map
  std::map<std::string, double> features;

  if (count == 0) {

    features["point"] = point_pto;
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


  // Compute both eigenvalues and eigenvectors
  if (SIZE > solver_thresh) {
    es.computeDirect(cov, true);
  } else {
    es.compute(cov, true);
  }


  Eigen::VectorXd eigenvalues = es.eigenvalues().real();

  // std::cout << "Eigenvalues:\n" << eigenvalues << std::endl;


  // Compute the verticality only if required
  double ver;
  if (Verticality) {

    Eigen::MatrixXd eigenvectors = es.eigenvectors().real();

    // std::cout << "Eigenvectors:\n" << eigenvectors << std::endl;


    // The normal vector is the eigenvector corresponding to the smallest eigenvalue
    // This represents the direction of least variation (surface normal)
    Eigen::Vector3d normal_vector = eigenvectors.col(0);

    // std::cout << "Normal vector:\n" << normal_vector << std::endl;

    // normal_vector.normalize();

    // std::cout << "Normalized vector:\n" << normal_vector << std::endl;

    // Define vertical reference direction
    Eigen::Vector3d vertical_reference(0, 0, 1);

    // Calculate verticality as alignment with vertical direction
    double dot_product = std::abs(normal_vector.dot(vertical_reference));

    // std::cout << "Normal vector:\n" << dot_product << std::endl;

    // For tree crowns: 1 = vertical surface, 0 = horizontal surface
    ver = 1.0 - dot_product;


  } else {
    ver = std::numeric_limits<double>::quiet_NaN();
  }

  //................................................................... add exceptions for eigenvalues which is: lambda1 >= lambda2 >= lambda3 >= 0.0
  // std::cout << point_pto << std::endl;

  double l1 = eigenvalues(2);
  double l2 = eigenvalues(1);
  double l3 = eigenvalues(0);

  // if (l1 < l2) {
  //   Rcpp::stop("We expect lambda1 to be equal or greater than lambda2!");
  // }
  // if (l2 < l3) {
  //   Rcpp::stop("We expect lambda2 to be equal or greater than lambda3!");
  // }
  // if (l3 < 0.0) {
  //   Rcpp::stop("We expect lambda3 to be greater than 0.0!");
  // }
  //...................................................................

  // Precompute values for reuse
  double lambda1 = (First_eigenvalue) ? l1 : std::numeric_limits<double>::quiet_NaN();
  double lambda2 = (Second_eigenvalue) ? l2 : std::numeric_limits<double>::quiet_NaN();
  double lambda3 = (Third_eigenvalue) ? l3 : std::numeric_limits<double>::quiet_NaN();
  double lambda123 = (Sum_of_eigenvalues) ? eigenvalues.sum() : std::numeric_limits<double>::quiet_NaN();

  features["point"] = point_pto;
  if (First_eigenvalue) features["First_eigenvalue"] = lambda1;
  if (Second_eigenvalue) features["Second_eigenvalue"] = lambda2;
  if (Third_eigenvalue) features["Third_eigenvalue"] = lambda3;
  if (Sum_of_eigenvalues) features["Sum_of_eigenvalues"] = lambda123;
  if (PCA_1) features["PCA1"] = (PCA_1) ? lambda1 / lambda123 : std::numeric_limits<double>::quiet_NaN();
  if (PCA_2) features["PCA2"] = (PCA_2) ? lambda2 / lambda123 : std::numeric_limits<double>::quiet_NaN();
  if (Anisotropy) features["Anisotropy"] = (Anisotropy) ? (lambda1 - lambda3) / lambda1 : std::numeric_limits<double>::quiet_NaN();
  if (Planarity) features["Planarity"] = (Planarity) ? (lambda2 - lambda3) / lambda1 : std::numeric_limits<double>::quiet_NaN();
  if (Linearity) features["Linearity"] = (Linearity) ? (lambda1 - lambda2) / lambda1 : std::numeric_limits<double>::quiet_NaN();
  if (Surface_variation) features["Surface_variation"] = (Surface_variation) ? lambda3 / lambda123 : std::numeric_limits<double>::quiet_NaN();
  if (Normal_change_rate) features["Normal_change_rate"] = (Normal_change_rate) ? lambda3 / lambda1 : std::numeric_limits<double>::quiet_NaN();
  if (Verticality) features["Verticality"] = ver;
  if (Number_of_points) features["Number_of_points"] = (Number_of_points) ? static_cast<double>(count) : std::numeric_limits<double>::quiet_NaN();
  if (omnivariance) features["omnivariance"] = (omnivariance) ? std::pow(lambda1 * lambda2 * lambda3, 1.0 / 3.0)  : std::numeric_limits<double>::quiet_NaN();
  if (eigenentropy) features["eigenentropy"] = (eigenentropy) ? -(lambda1 * std::log(lambda1) + lambda2 * std::log(lambda2) + lambda3 * std::log(lambda3)) : std::numeric_limits<double>::quiet_NaN();
  if (surface_density) features["surface_density"] = (surface_density) ? static_cast<double>(count) / (4.0 * M_PI * std::pow(dist, 2.0))  : std::numeric_limits<double>::quiet_NaN();
  if (volume_density) features["volume_density"] = (volume_density) ? static_cast<double>(count) / ((4.0/3.0) * M_PI * std::pow(dist, 3.0))  : std::numeric_limits<double>::quiet_NaN();

  return features;
}





// Rcpp function that takes a matrix of 4 columns (point, x, y, z) as input and a list of lists, where each sublist is a vector of indices.
// Then the function iterates over each sublist and subsets the matrix by each vector of indices.
// We account also for the fact that the indexing in the vectors start from 1 whereas in Rcpp from 0

// [[Rcpp::export]]
Rcpp::List subset_matrix_by_indices(const Eigen::MatrixXd& mat,
                                    Rcpp::List index_lists,
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
                                    int solver_thresh = 50000) {

  // Check if matrix has 4 columns
  if (mat.cols() != 4) {
    Rcpp::stop("Matrix must have exactly 4 columns");
  }

  int n_sublists = index_lists.size();
  Rcpp::List result(n_sublists);

  // Iterate over each sublist of indices
  for (int i = 0; i < n_sublists; i++) {

    // Extract the current vector of indices
    Rcpp::IntegerVector indices = index_lists[i];
    int n_indices = indices.size();

    double point_i = mat(i, 0);

    // Handle empty index vector - add NULL to result
    if (n_indices == 0) {

      std::map<std::string, double> features_null;

      features_null["point"] = point_i;

      if (First_eigenvalue) features_null["First_eigenvalue"] = std::numeric_limits<double>::quiet_NaN();
      if (Second_eigenvalue) features_null["Second_eigenvalue"] = std::numeric_limits<double>::quiet_NaN();
      if (Third_eigenvalue) features_null["Third_eigenvalue"] = std::numeric_limits<double>::quiet_NaN();
      if (Sum_of_eigenvalues) features_null["Sum_of_eigenvalues"] = std::numeric_limits<double>::quiet_NaN();
      if (PCA_1) features_null["PCA1"] = std::numeric_limits<double>::quiet_NaN();
      if (PCA_2) features_null["PCA2"] = std::numeric_limits<double>::quiet_NaN();
      if (Anisotropy) features_null["Anisotropy"] = std::numeric_limits<double>::quiet_NaN();
      if (Planarity) features_null["Planarity"] = std::numeric_limits<double>::quiet_NaN();
      if (Linearity) features_null["Linearity"] = std::numeric_limits<double>::quiet_NaN();
      if (Surface_variation) features_null["Surface_variation"] = std::numeric_limits<double>::quiet_NaN();
      if (Normal_change_rate) features_null["Normal_change_rate"] = std::numeric_limits<double>::quiet_NaN();
      if (Verticality) features_null["Verticality"] = std::numeric_limits<double>::quiet_NaN();
      if (Number_of_points) features_null["Number_of_points"] = std::numeric_limits<double>::quiet_NaN();
      if (omnivariance) features_null["omnivariance"] = std::numeric_limits<double>::quiet_NaN();
      if (eigenentropy) features_null["eigenentropy"] = std::numeric_limits<double>::quiet_NaN();
      if (surface_density) features_null["surface_density"] = std::numeric_limits<double>::quiet_NaN();
      if (volume_density) features_null["volume_density"] = std::numeric_limits<double>::quiet_NaN();
      result[i] = features_null;

      continue;
    }

    // Check for valid indices (must be >= 1 and <= number of rows)
    for (int j = 0; j < n_indices; j++) {
      if (indices[j] < 1 || indices[j] > mat.rows()) {
        Rcpp::stop("Index out of bounds: indices must be between 1 and %d", mat.rows());
      }
    }

    // Create subset matrix
    Eigen::MatrixXd subset(n_indices, 4);

    // Fill the subset matrix
    for (int row = 0; row < n_indices; row++) {
      int original_row = indices[row] - 1;                // Convert from 1-based to 0-based indexing
      for (int col = 0; col < 4; col++) {
        subset(row, col) = mat(original_row, col);
      }
    }

    // Extract columns as mapped vectors from subset matrix
    const Eigen::Map<Eigen::VectorXd> point(subset.col(0).data(), n_indices);
    const Eigen::Map<Eigen::VectorXd> x(subset.col(1).data(), n_indices);
    const Eigen::Map<Eigen::VectorXd> y(subset.col(2).data(), n_indices);
    const Eigen::Map<Eigen::VectorXd> z(subset.col(3).data(), n_indices);

    // Extract values from row i of the original matrix (assuming i is the current sublist index)
    double x_pto = mat(i, 1);
    double y_pto = mat(i, 2);
    double z_pto = mat(i, 3);

    std::map<std::string, double> res_geom = geometric_features_point(x,
                                                                      y,
                                                                      z,
                                                                      x_pto,
                                                                      y_pto,
                                                                      z_pto,
                                                                      point_i,
                                                                      dist,
                                                                      First_eigenvalue,
                                                                      Second_eigenvalue,
                                                                      Third_eigenvalue,
                                                                      Sum_of_eigenvalues,
                                                                      PCA_1,
                                                                      PCA_2,
                                                                      Anisotropy,
                                                                      Planarity,
                                                                      Linearity,
                                                                      Surface_variation,
                                                                      Normal_change_rate,
                                                                      Verticality,
                                                                      Number_of_points,
                                                                      omnivariance,
                                                                      eigenentropy,
                                                                      surface_density,
                                                                      volume_density,
                                                                      solver_thresh);



    result[i] = res_geom;
  }

  return result;
}
