// #define STRICT_R_HEADERS
// #include <Rcpp.h>
// #include <RcppEigen.h>
//
// // [[Rcpp::depends(RcppEigen)]]
//
// using namespace Rcpp;
// using namespace std;
//
//
// using Eigen::Map;
// using Eigen::MatrixXd;
// using Eigen::VectorXd;
//
//
// double  getEuclideanDistance(double x1, double y1, double z1, double x2, double y2, double z2) {
//
//   double dx, dy, dz, d;
//
//   dx = x2 - x1;
//   dy = y2 - y1;
//   dz = z2 - z1;
//
//   d = dx * dx + dy * dy + dz * dz;
//
//   return sqrt(d);
//
// }
//
//
// // double computeLastEigenvalue(Eigen::MatrixXd matrix) {
// //   Eigen::EigenSolver<Eigen::MatrixXd> solver(matrix);
// //   Eigen::VectorXd eigenvalues = solver.eigenvalues().real();
// //   return eigenvalues(eigenvalues.size() - 1);
// // }
//
//
//
// double ncr_double(Eigen::VectorXd x, Eigen::VectorXd y, Eigen::VectorXd z, int pto, double dist) {
//
//   double d = 0.0;
//   double dd = 0.0;
//
//   vector<double> xx(0);
//   vector<double> yy(0);
//   vector<double> zz(0);
//
//   for (int i=0; i<x.size(); i++){
//
//     d = getEuclideanDistance(x[pto], y[pto], z[pto], x[i], y[i], z[i]);
//
//     if (d <= dist){
//
//       xx.push_back(x[i]);
//       yy.push_back(y[i]);
//       zz.push_back(z[i]);
//
//     } else {
//
//       continue;
//
//     }
//
//   }
//
//   Eigen::MatrixXd m = Eigen::MatrixXd::Zero(xx.size(), 3);
//
//   for (int j=0; j<xx.size(); j++){
//
//     m(j,0) = xx[j];
//     m(j,1) = yy[j];
//     m(j,2) = zz[j];
//
//   }
//
//   Eigen::MatrixXd centered = m.rowwise() - m.colwise().mean();
//   Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(m.rows() - 1);
//
//   Eigen::SelfAdjointEigenSolver<MatrixXd> es;
//   es.compute(cov, /* computeEigenvectors = */ false);
//
//   Eigen::VectorXd kk = es.eigenvalues().real();
//
//   dd = kk(0, 0) / kk.sum();
//   //dd = round(dd * 10000.0) / 10000.0;
//
//   //std::cout << "Here is the matrix" << m << std::endl;
//
//   return dd;
//
// }
//
//
// double ver_double(Eigen::VectorXd x, Eigen::VectorXd y, Eigen::VectorXd z, int pto, double dist) {
//
//   double d = 0.0;
//   double dd = 0.0;
//   double lastEigenValue = 0.0;
//
//   vector<double> xx(0);
//   vector<double> yy(0);
//   vector<double> zz(0);
//
//   for (int i=0; i<x.size(); i++){
//
//     d = getEuclideanDistance(x[pto], y[pto], z[pto], x[i], y[i], z[i]);
//
//     if (d <= dist){
//
//       xx.push_back(x[i]);
//       yy.push_back(y[i]);
//       zz.push_back(z[i]);
//
//     } else {
//
//       continue;
//
//     }
//
//   }
//
//   Eigen::MatrixXd m = Eigen::MatrixXd::Zero(xx.size(), 3);
//
//   for (int j=0; j<xx.size(); j++){
//
//     m(j,0) = xx[j];
//     m(j,1) = yy[j];
//     m(j,2) = zz[j];
//
//   }
//
//   Eigen::MatrixXd centered = m.rowwise() - m.colwise().mean();
//   Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(m.rows() - 1);
//
//   Eigen::SelfAdjointEigenSolver<MatrixXd> es;
//   es.compute(cov, /* computeEigenvectors = */ false);
//
//   Eigen::VectorXd eigenvalues = es.eigenvalues().real();
//   Eigen::MatrixXd eigenvectors = es.eigenvectors().real();
//
//   // Find the eigenvector corresponding to the smallest eigenvalue
//   int minEigenIndex = 0;
//   for (int i = 1; i < 3; ++i) {
//     if (eigenvalues(i) < eigenvalues(minEigenIndex)) {
//       minEigenIndex = i;
//     }
//   }
//
//   Eigen::Vector3d normalVector = eigenvectors.col(minEigenIndex);
//
//   // Normalize the normal vector to obtain a unit vector
//   normalVector.normalize();
//
//
//   // dd = 1 - lastEigenValue / kk.sum();
//   dd = 1 - abs(normalVector[2]);
//
//   // std::cout << "Here are the eigenvalues:\n" << kk << std::endl;
//   // std::cout << "Here is the last eigenvalue:\n" << lastEigenValue << std::endl;
//
//   return dd;
//
// }
//
//
//
// std::vector<double> geometric_features_point(Eigen::VectorXd x, Eigen::VectorXd y, Eigen::VectorXd z, int pto, double dist) {
//
//   double d = 0.0;
//
//   std::vector<double> result;
//
//   vector<double> xx(0);
//   vector<double> yy(0);
//   vector<double> zz(0);
//
//   for (int i=0; i<x.size(); i++){
//
//     d = getEuclideanDistance(x[pto], y[pto], z[pto], x[i], y[i], z[i]);
//
//     if (d <= dist){
//
//       xx.push_back(x[i]);
//       yy.push_back(y[i]);
//       zz.push_back(z[i]);
//
//     } else {
//
//       continue;
//
//     }
//
//   }
//
//   Eigen::MatrixXd m = Eigen::MatrixXd::Zero(xx.size(), 3);
//
//   for (int j=0; j<xx.size(); j++){
//
//     m(j,0) = xx[j];
//     m(j,1) = yy[j];
//     m(j,2) = zz[j];
//
//   }
//
//   Eigen::MatrixXd centered = m.rowwise() - m.colwise().mean();
//   Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(m.rows() - 1);
//
//   Eigen::SelfAdjointEigenSolver<MatrixXd> es;
//   es.compute(cov, /* computeEigenvectors = */ false);
//
//   Eigen::VectorXd eigenvalues = es.eigenvalues().real();
//   Eigen::MatrixXd eigenvectors = es.eigenvectors().real();
//
//   // Find the eigenvector corresponding to the smallest eigenvalue
//   int minEigenIndex = 0;
//   for (int i = 1; i < 3; ++i) {
//     if (eigenvalues(i) < eigenvalues(minEigenIndex)) {
//       minEigenIndex = i;
//     }
//   }
//
//   Eigen::Vector3d normalVector = eigenvectors.col(minEigenIndex);
//
//   // Normalize the normal vector to obtain a unit vector
//   normalVector.normalize();
//
//   // Geometric features:
//
//   double lambda1 = eigenvalues(2); // First eigenvalue
//   double lambda2 = eigenvalues(1); // Second eigenvalue
//   double lambda3 = eigenvalues(0); // Third eigenvalue
//   double lambda123 = eigenvalues.sum(); // Sum of eigenvalues
//   double PCA1 = lambda1 / eigenvalues.sum(); // PCA1
//   double PCA2 = lambda2 / eigenvalues.sum(); // PCA2
//   double ani = (lambda1 - lambda3) / lambda1; // Anisotropy
//   double pla = (lambda2 - lambda3) / lambda1; // Planarity
//   double lin = (lambda1 - lambda2) / lambda1; // Linearity
//   double sur = lambda3 / eigenvalues.sum(); // Surface variation
//   double sph = lambda3 / lambda1; // Normal change rate
//   // double ver = 1 - (abs(normalVector[2]) / normalVector.sum()); // Verticality
//   double ver = 1 - abs(normalVector[2]); // Verticality
//   double pts = m.rows(); // Normal change rate
//
//   // std::cout << "Here is the matrix" << kk << std::endl;
//
//   result.push_back(lambda1);
//   result.push_back(lambda2);
//   result.push_back(lambda3);
//   result.push_back(lambda123);
//   result.push_back(PCA1);
//   result.push_back(PCA2);
//   result.push_back(ani);
//   result.push_back(pla);
//   result.push_back(lin);
//   result.push_back(sur);
//   result.push_back(sph);
//   result.push_back(ver);
//   result.push_back(pts);
//
//   return result;
//
// }
//
//
// // [[Rcpp::export]]
// DataFrame ver_point_cloud_double(const Eigen::MatrixXd & m) {
//
//   Eigen::VectorXd point = m.col(0);
//   Eigen::VectorXd x = m.col(1);
//   Eigen::VectorXd y = m.col(2);
//   Eigen::VectorXd z = m.col(3);
//
//   double firstX = x.minCoeff();
//   double lastX = x.maxCoeff();
//
//   double firstY = y.minCoeff();
//   double lastY = y.maxCoeff();
//
//   double firstZ = z.minCoeff();
//   double lastZ = z.maxCoeff();
//
//   int p = 0;
//   double k = 0.0;
//   double dist = 0.05;
//
//   vector<double> kk(0);
//   vector<int> kkk(0);
//
//   for(int q = 0; q < x.size(); q++) {
//
//     if(x[q] > lastX - dist || x[q] < firstX + dist || y[q] > lastY - dist || y[q] < firstY + dist || z[q] > lastZ - dist || z[q] < firstZ + dist){
//
//       continue;
//
//     } else {
//
//       k = ver_double(x, y, z, q, dist);
//       p = point[q];
//
//     } if(isnan(k)){
//
//       continue;
//
//     } else {
//
//       kk.push_back(k);
//       kkk.push_back(p);
//
//     }
//
//   }
//
//   if(kk.size() == 0){
//
//     kk.push_back(9999);
//     kkk.push_back(p);
//
//   }
//
//   return DataFrame::create(Named("point")=kkk,
//                            Named("ver")=kk);
//
// }
//
//
// // [[Rcpp::export]]
// DataFrame ncr_point_cloud_double(const Eigen::MatrixXd & m) {
//
//   Eigen::VectorXd point = m.col(0);
//   Eigen::VectorXd x = m.col(1);
//   Eigen::VectorXd y = m.col(2);
//   Eigen::VectorXd z = m.col(3);
//
//   double firstX = x.minCoeff();
//   double lastX = x.maxCoeff();
//
//   double firstY = y.minCoeff();
//   double lastY = y.maxCoeff();
//
//   double firstZ = z.minCoeff();
//   double lastZ = z.maxCoeff();
//
//   int p = 0;
//   double k = 0.0;
//   double dist = 0.05;
//
//   vector<double> kk(0);
//   vector<int> kkk(0);
//
//   for(int q = 0; q < x.size(); q++) {
//
//     if(x[q] > lastX - dist || x[q] < firstX + dist || y[q] > lastY - dist || y[q] < firstY + dist || z[q] > lastZ - dist || z[q] < firstZ + dist){
//
//       continue;
//
//     } else {
//
//       k = ncr_double(x, y, z, q, dist);
//       p = point[q];
//
//     } if(isnan(k)){
//
//       continue;
//
//     } else {
//
//       kk.push_back(k);
//       kkk.push_back(p);
//
//     }
//
//   }
//
//   if(kk.size() == 0){
//
//      kk.push_back(9999);
//      kkk.push_back(p);
//
//   }
//
//   return DataFrame::create(Named("point")=kkk,
//                            Named("ncr")=kk);
//
// }
//
//
// // [[Rcpp::export]]
// DataFrame geometric_features(const Eigen::MatrixXd & m, double dist) {
//
//   Eigen::VectorXd point = m.col(0);
//   Eigen::VectorXd x = m.col(1);
//   Eigen::VectorXd y = m.col(2);
//   Eigen::VectorXd z = m.col(3);
//
//   double firstX = x.minCoeff();
//   double lastX = x.maxCoeff();
//
//   double firstY = y.minCoeff();
//   double lastY = y.maxCoeff();
//
//   double firstZ = z.minCoeff();
//   double lastZ = z.maxCoeff();
//
//   // double dist = 0.05;
//
//   vector<double> lambda1(0);
//   vector<double> lambda2(0);
//   vector<double> lambda3(0);
//   vector<double> lambda123(0);
//   vector<double> PCA1(0);
//   vector<double> PCA2(0);
//   vector<double> ani(0);
//   vector<double> pla(0);
//   vector<double> lin(0);
//   vector<double> sur(0);
//   vector<double> sph(0);
//   vector<double> ver(0);
//   vector<double> pts(0);
//
//   vector<double> k(0);
//   vector<int> kk(0);
//
//   for(int q = 0; q < x.size(); q++) {
//
//     if(x[q] > lastX - dist || x[q] < firstX + dist || y[q] > lastY - dist || y[q] < firstY + dist || z[q] > lastZ - dist || z[q] < firstZ + dist){
//
//       continue;
//
//     }
//
//     k = geometric_features_point(x, y, z, q, dist);
//
//     bool hasNaN = false;
//     for (double val : k) {
//       if (std::isnan(val)) {
//         hasNaN = true;
//         break;
//       }
//     }
//
//     if (hasNaN) {
//       continue;
//     }
//
//       lambda1.push_back(k[0]);
//       lambda2.push_back(k[1]);
//       lambda3.push_back(k[2]);
//       lambda123.push_back(k[3]);
//       PCA1.push_back(k[4]);
//       PCA2.push_back(k[5]);
//       ani.push_back(k[6]);
//       pla.push_back(k[7]);
//       lin.push_back(k[8]);
//       sur.push_back(k[9]);
//       sph.push_back(k[10]);
//       ver.push_back(k[11]);
//       pts.push_back(k[12]);
//       kk.push_back(point[q]);
//
//     }
//
//   if(lambda1.empty()){
//      lambda1.push_back(9999);
//      kk.push_back(point[0]);
//   }
//
//   if(lambda2.empty()){
//      lambda2.push_back(9999);
//   }
//
//   if(lambda3.empty()){
//      lambda3.push_back(9999);
//   }
//
//   if(lambda123.empty()){
//      lambda123.push_back(9999);
//   }
//
//   if(PCA1.empty()){
//     PCA1.push_back(9999);
//   }
//
//   if(PCA2.empty()){
//     PCA2.push_back(9999);
//   }
//
//   if(ani.empty()){
//      ani.push_back(9999);
//   }
//
//   if(pla.empty()){
//      pla.push_back(9999);
//   }
//
//   if(lin.empty()){
//      lin.push_back(9999);
//   }
//
//   if(sur.empty()){
//      sur.push_back(9999);
//   }
//
//   if(sph.empty()){
//      sph.push_back(9999);
//   }
//
//   if(ver.empty()){
//      ver.push_back(9999);
//   }
//
//   if(pts.empty()){
//      pts.push_back(9999);
//   }
//
//   // Create and return DataFrame
//   return DataFrame::create(Named("point")=kk,
//                            Named("first_eigenvalue")=lambda1,
//                            Named("second_eigenvalue")=lambda2,
//                            Named("third_eigenvalue")=lambda3,
//                            Named("sum_eigenvalues")=lambda123,
//                            Named("PCA1")=PCA1,
//                            Named("PCA2")=PCA2,
//                            Named("anisotropy")=ani,
//                            Named("planarity")=pla,
//                            Named("linearity")=lin,
//                            Named("surface_variation")=sur,
//                            Named("sphericity")=sph,
//                            Named("verticality")=ver,
//                            Named("number_neighbors")=pts);
//
// }
