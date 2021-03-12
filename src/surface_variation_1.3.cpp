
#include <Rcpp.h>
#include <RcppEigen.h>


// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;


using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;


double  getEuclideanDistance(double x1, double y1, double z1, double x2, double y2, double z2) {

  double dx, dy, dz, d;

  dx = x2 - x1;
  dy = y2 - y1;
  dz = z2 - z1;

  d = dx * dx + dy * dy + dz * dz;

  return sqrt(d);

}


double ncr_double(Eigen::VectorXd x, Eigen::VectorXd y, Eigen::VectorXd z, int pto, double dist) {

  double d = 0.0;
  double dd = 0.0;

  vector<double> xx(0);
  vector<double> yy(0);
  vector<double> zz(0);

  for (int i=0; i<x.size(); i++){

    d = getEuclideanDistance(x[pto], y[pto], z[pto], x[i], y[i], z[i]);

    if (d <= dist){

      xx.push_back(x[i]);
      yy.push_back(y[i]);
      zz.push_back(z[i]);

    } else {

      continue;

    }

  }

  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(xx.size(), 3);

  for (int j=0; j<xx.size(); j++){

    m(j,0) = xx[j];
    m(j,1) = yy[j];
    m(j,2) = zz[j];

  }

  Eigen::MatrixXd centered = m.rowwise() - m.colwise().mean();
  Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(m.rows() - 1);

  Eigen::SelfAdjointEigenSolver<MatrixXd> es;
  es.compute(cov, /* computeEigenvectors = */ false);

  Eigen::VectorXd kk = es.eigenvalues().real();

  dd = kk(0, 0) / kk.sum();
  //dd = round(dd * 10000.0) / 10000.0;

  //std::cout << "Here is the matrix" << m << std::endl;

  return dd;

}

// [[Rcpp::export]]
DataFrame ncr_point_cloud_double(const Eigen::MatrixXd & m) {

  Eigen::VectorXd point = m.col(0);
  Eigen::VectorXd x = m.col(1);
  Eigen::VectorXd y = m.col(2);
  Eigen::VectorXd z = m.col(3);

  double firstX = x.minCoeff();
  double lastX = x.maxCoeff();

  double firstY = y.minCoeff();
  double lastY = y.maxCoeff();

  double firstZ = z.minCoeff();
  double lastZ = z.maxCoeff();

  int p = 0;
  double k = 0.0;
  double dist = 0.05;

  vector<double> kk(0);
  vector<int> kkk(0);

  for(int q = 0; q < x.size(); q++) {

    if(x[q] > lastX - dist || x[q] < firstX + dist || y[q] > lastY - dist || y[q] < firstY + dist || z[q] > lastZ - dist || z[q] < firstZ + dist){

      continue;

    } else {

      k = ncr_double(x, y, z, q, dist);
      p = point[q];

    } if(isnan(k)){

      continue;

    } else {

      kk.push_back(k);
      kkk.push_back(p);

    }

  }

  if(kk.size() == 0){

     kk.push_back(9999);
     kkk.push_back(p);

  }

  return DataFrame::create(Named("point")=kkk,
                           Named("ncr")=kk);

}
