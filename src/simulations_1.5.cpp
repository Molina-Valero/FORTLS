#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double weighted_mean_arit(std::vector<double> x, std::vector<double> w) {
  // Initialize these to zero
  double total_w = 0;
  double total_xw = 0;

  // Set n to the size of x
  int n = x.size();

  // Specify the for loop arguments
  for(int i = 0; i < n; i++) {
    // Add ith weight
    total_w += w[i];

    // Add the ith data value times the ith weight
    total_xw += x[i] * w[i];

  }

  // Return the total product divided by the total weight
  return total_xw / total_w;
}

// [[Rcpp::export]]
double weighted_mean_sqrt(std::vector<double> x, std::vector<double> w) {
  // Initialize these to zero
  double total_w = 0;
  double total_xw = 0;

  // Set n to the size of x
  int n = x.size();

  // Specify the for loop arguments
  for(int i = 0; i < n; i++) {
    // Add ith weight
    total_w += w[i];
    // Add the ith data value times the ith weight
    total_xw += pow(x[i], 2.0) * w[i];
  }

  // Return the total product divided by the total weight
  return sqrt(total_xw / total_w);
}

// [[Rcpp::export]]
double weighted_mean_geom(std::vector<double> x, std::vector<double> w) {
  // Initialize these to zero
  double total_w = 0;
  double total_xw = 1;
  double exp = 0;

  // Set n to the size of x
  int n = x.size();

  // Specify the for loop arguments
  for(int i = 0; i < n; i++) {
    //Add ith weight
    total_w += w[i];

    // Add the ith data value times the ith weight
    total_xw *= pow(x[i], w[i]);

    exp = 1 / total_w;

  }

  // Return the total product divided by the total weight
  return pow(total_xw, exp);
}

// [[Rcpp::export]]
double weighted_mean_harm(std::vector<double> x, std::vector<double> w) {
  // Initialize these to zero
  double total_w = 0;
  double total_xw = 0;
  double y = 0;

  // Set n to the size of x
  int n = x.size();

  // Specify the for loop arguments
  for(int i = 0; i < n; i++) {
    // Add ith weight
    total_w += w[i];
    // Add the ith data value times the ith weight

    if(w[i] == 0.0){

      y = 0.0;

    } else {

      y = w[i] / x[i];

    }

    total_xw += y;

  }

  // Return the total product divided by the total weight
  return total_w / total_xw;
}



template<typename T>
static inline double Lerp(T v0, T v1, T t)
{
  return (1 - t)*v0 + t*v1;
}

template<typename T>
static inline std::vector<T> Quantile(const std::vector<T>& inData,
                                      const std::vector<T>& probs)
{
  if (inData.empty())
  {
    return std::vector<T>();
  }

  if (1 == inData.size())
  {
    return std::vector<T>(1, inData[0]);
  }

  std::vector<T> data = inData;
  std::sort(data.begin(), data.end());
  std::vector<T> quantiles;

  for (size_t i = 0; i < probs.size(); ++i)
  {
    T poi = Lerp<T>(-0.5, data.size() - 0.5, probs[i]);

    size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
    size_t right = std::min(int64_t(std::ceil(poi)),
                            int64_t(data.size() - 1));

    T datLeft = data.at(left);
    T datRight = data.at(right);

    T quantile = Lerp<T>(datLeft, datRight, poi - left);

    quantiles.push_back(quantile);
  }

  return quantiles;
}


// [[Rcpp::export]]
DataFrame height_perc_cpp(std::vector<double> rho_seq, std::vector<double> z,
                          std::vector<double> rho) {

  int n = rho_seq.size();

  NumericVector P01(n);
  NumericVector P05(n);
  NumericVector P10(n);
  NumericVector P20(n);
  NumericVector P25(n);
  NumericVector P30(n);
  NumericVector P40(n);
  NumericVector P50(n);
  NumericVector P60(n);
  NumericVector P70(n);
  NumericVector P75(n);
  NumericVector P80(n);
  NumericVector P90(n);
  NumericVector P95(n);
  NumericVector P99(n);
  NumericVector P999(n);

  vector<double> quartiles(14);

  int n_z = z.size();
  vector<double> percentile(n_z);

  vector<double>::iterator pos;

  for(int i = 0; i < n; i++){

    pos = find_if(rho.begin(), rho.end(),
                  std::bind(std::greater<double>(), std::placeholders::_1, rho_seq[i]));

    percentile.assign(z.begin(), z.begin() + distance(rho.begin(), pos));

    quartiles = Quantile<double>(percentile, {0.01, 0.05, 0.10,
                                              0.20, 0.25, 0.30,
                                              0.40, 0.50, 0.60,
                                              0.70, 0.75, 0.80,
                                              0.90, 0.95, 0.99, 0.999});

    P01[i] = quartiles[0];
    P05[i] = quartiles[1];
    P10[i] = quartiles[2];
    P20[i] = quartiles[3];
    P25[i] = quartiles[4];
    P30[i] = quartiles[5];
    P40[i] = quartiles[6];
    P50[i] = quartiles[7];
    P60[i] = quartiles[8];
    P70[i] = quartiles[9];
    P75[i] = quartiles[10];
    P80[i] = quartiles[11];
    P90[i] = quartiles[12];
    P95[i] = quartiles[13];
    P99[i] = quartiles[14];
    P999[i] = quartiles[15];

  }

  return DataFrame::create(Named("rho_seq") = rho_seq,
                           Named("P01") = P01, Named("P05") = P05,
                           Named("P10") = P10, Named("P20") = P20,
                           Named("P25") = P25, Named("P30") = P30,
                           Named("P40") = P40, Named("P50") = P50,
                           Named("P60") = P60, Named("P70") = P70,
                           Named("P75") = P75, Named("P80") = P80,
                           Named("P90") = P90, Named("P95") = P95,
                           Named("P99") = P99, Named("P99.9") = P999);

}

vector<double> order_cpp(vector<double> x, vector<double> y) {

  std::vector<double> d(x.size());
  std::vector<double> h(y.size());

  // Declaring vector of pairs

  // create a vector of length of the smaller vector
  std::vector<std::pair<double, double>> target( x.size() < y.size() ? x.size() : y.size() );

  // Entering values in vector of pairs
  for (unsigned i = 0; i < target.size(); i++)
    target[i] = std::make_pair(x[i], y[i]);


  // Using modified sort() function to sort
  sort(target.rbegin(), target.rend());

  for (int i = 0; i < target.size(); i++){

    d[i] = target[i].first;
    h[i] = target[i].second;

  }

  return h;

}

// [[Rcpp::export]]
DataFrame fixed_area_cpp(std::vector<double> radius_seq,
                       std::vector<double> hdist, std::vector<double> d,
                       std::vector<double> h, double num) {

  int n = radius_seq.size();

  NumericVector d_dom_arit(n);
  NumericVector d_dom_sqrt(n);
  NumericVector d_dom_geom(n);
  NumericVector d_dom_harm(n);

  NumericVector h_dom_arit(n);
  NumericVector h_dom_sqrt(n);
  NumericVector h_dom_geom(n);
  NumericVector h_dom_harm(n);

  int n_d = d.size();

  vector<double> x(n_d);
  vector<double> y(n_d);
  vector<double> h_0(n_d);

  vector<double> weights(n_d, 1);
  vector<double> factor(n_d);

  vector<double>::iterator pos;

  double ef = 1;

  for(int i = 0; i < n; i++){

    pos = find_if(hdist.begin(), hdist.end(),
                  std::bind(std::greater<double>(), std::placeholders::_1, radius_seq[i]));

    x.assign(d.begin(), d.begin() + distance(hdist.begin(), pos));
    y.assign(h.begin(), h.begin() + distance(hdist.begin(), pos));

    h_0 = order_cpp(x, y);

    sort(x.begin(), x.end(), greater<double>());

    factor.assign(weights.begin(), weights.begin()
                                   + distance(hdist.begin(), pos));
    ef = 10000 / (M_PI * pow(radius_seq[i], 2.0));
    transform(factor.begin(), factor.end(), factor.begin(),
              std::bind(multiplies<double>(), std::placeholders::_1, ef));
    partial_sum(factor.begin(), factor.end(), factor.begin());

    vector<double> factor_1(factor);
    factor_1.insert(factor_1.begin(), 0.0);

    for(int j = 0; j < factor.size(); j++){

      if(factor[j] <= num) {

        factor[j] = 1.0;

      } else if ((factor[j] > num) & (factor_1[j] < num)) {

        factor[j] = (num - factor_1[j]) / num;

      } else {

        factor[j] = 0.0;

      }

    }

    d_dom_arit[i] = weighted_mean_arit(x, factor);
    d_dom_sqrt[i] = weighted_mean_sqrt(x, factor);
    d_dom_geom[i] = weighted_mean_geom(x, factor);
    d_dom_harm[i] = weighted_mean_harm(x, factor);

    h_dom_arit[i] = weighted_mean_arit(h_0, factor);
    h_dom_sqrt[i] = weighted_mean_sqrt(h_0, factor);
    h_dom_geom[i] = weighted_mean_geom(h_0, factor);
    h_dom_harm[i] = weighted_mean_harm(h_0, factor);

  }

  return DataFrame::create(Named("radius") = radius_seq,
                           Named("d.0") = d_dom_arit,
                           Named("dg.0") = d_dom_sqrt,
                           Named("dgeom.0") = d_dom_geom,
                           Named("dharm.0") = d_dom_harm,
                           Named("h.0") = h_dom_arit,
                           Named("hg.0") = h_dom_sqrt,
                           Named("hgeom.0") = h_dom_geom,
                           Named("hharm.0") = h_dom_harm);

}

// [[Rcpp::export]]
DataFrame k_tree_cpp(std::vector<double> k_seq,
                     std::vector<double> radius_seq, std::vector<double> k,
                     std::vector<double> d, std::vector<double> h,
                     double num) {

  int n = k_seq.size();

  NumericVector d_dom_arit(n);
  NumericVector d_dom_sqrt(n);
  NumericVector d_dom_geom(n);
  NumericVector d_dom_harm(n);

  NumericVector h_dom_arit(n);
  NumericVector h_dom_sqrt(n);
  NumericVector h_dom_geom(n);
  NumericVector h_dom_harm(n);

  int n_d = d.size();

  vector<double> x(n_d);
  vector<double> y(n_d);
  vector<double> h_0(n_d);

  vector<double> weights(n_d, 1);
  vector<double> factor(n_d);

  vector<double>::iterator pos;

  double ef = 1;

  for(int i = 0; i < n; i++){

    pos = find_if(k.begin(), k.end(),
                  std::bind(std::greater<double>(), std::placeholders::_1, k_seq[i]));


    x.assign(d.begin(), d.begin() + distance(k.begin(), pos));
    y.assign(h.begin(), h.begin() + distance(k.begin(), pos));

    h_0 = order_cpp(x, y);

    sort(x.begin(), x.end(), greater<double>());

    factor.assign(weights.begin(), weights.begin()
                                   + distance(k.begin(), pos));
    ef = 10000 / (M_PI * pow(radius_seq[i], 2.0));
    transform(factor.begin(), factor.end(), factor.begin(),
              std::bind(multiplies<double>(), std::placeholders::_1, ef));
    partial_sum(factor.begin(), factor.end(), factor.begin());

    vector<double> factor_1(factor);
    factor_1.insert(factor_1.begin(), 0.0);

    for(int j = 0; j < factor.size(); j++){

      if(factor[j] <= num) {

        factor[j] = 1.0;

      } else if ((factor[j] > num) & (factor_1[j] < num)) {

        factor[j] = (num - factor_1[j]) / num;

      } else {

        factor[j] = 0.0;

      }

    }

    d_dom_arit[i] = weighted_mean_arit(x, factor);
    d_dom_sqrt[i] = weighted_mean_sqrt(x, factor);
    d_dom_geom[i] = weighted_mean_geom(x, factor);
    d_dom_harm[i] = weighted_mean_harm(x, factor);

    h_dom_arit[i] = weighted_mean_arit(h_0, factor);
    h_dom_sqrt[i] = weighted_mean_sqrt(h_0, factor);
    h_dom_geom[i] = weighted_mean_geom(h_0, factor);
    h_dom_harm[i] = weighted_mean_harm(h_0, factor);

  }

  return DataFrame::create(Named("k") = k_seq,
                           Named("d.0") = d_dom_arit,
                           Named("dg.0") = d_dom_sqrt,
                           Named("dgeom.0") = d_dom_geom,
                           Named("dharm.0") = d_dom_harm,
                           Named("h.0") = h_dom_arit,
                           Named("hg.0") = h_dom_sqrt,
                           Named("hgeom.0") = h_dom_geom,
                           Named("hharm.0") = h_dom_harm);

}


// [[Rcpp::export]]
DataFrame angle_count_cpp(std::vector<double> baf_seq,
                          std::vector<double> baf, std::vector<double> d,
                          std::vector<double> h, double num) {

  int n = baf_seq.size();

  NumericVector d_dom_arit(n);
  NumericVector d_dom_sqrt(n);
  NumericVector d_dom_geom(n);
  NumericVector d_dom_harm(n);

  NumericVector h_dom_arit(n);
  NumericVector h_dom_sqrt(n);
  NumericVector h_dom_geom(n);
  NumericVector h_dom_harm(n);

  int n_d = d.size();

  vector<double> x(n_d);
  vector<double> y(n_d);
  vector<double> h_0(n_d);

  vector<double> weights(n_d, 1);
  vector<double> factor(n_d);

  vector<double>::iterator pos;

  double ef = 1;


  for(int i = 0; i < n; i++){

    pos = find_if(baf.begin(), baf.end(),
                  std::bind(std::less<double>(), std::placeholders::_1, baf_seq[i]));

    x.assign(d.begin(), d.begin() + distance(baf.begin(), pos));
    y.assign(h.begin(), h.begin() + distance(baf.begin(), pos));

    h_0 = order_cpp(x, y);

    sort(x.begin(), x.end(), greater<double>());

    factor.assign(weights.begin(), weights.begin()
                                   + distance(baf.begin(), pos));
    ef = baf_seq[i];
    transform(factor.begin(), factor.end(), x.begin(), factor.begin(),
              [ef](double factor, double x) { return ef / ((M_PI / 4) *
                                              pow(x, 2.0)); });

    for(int q = 0; q < factor.size(); q++){

      if (isinf(factor[q])){

        factor[q] = 0.0;

      }

    }

    partial_sum(factor.begin(), factor.end(), factor.begin());

    vector<double> factor_1(factor);
    factor_1.insert(factor_1.begin(), 0.0);

    for(int j = 0; j < factor.size(); j++){

      if(factor[j] <= num) {

        factor[j] = 1.0;

      } else if ((factor[j] > num) & (factor_1[j] < num)) {

        factor[j] = (num - factor_1[j]) / num;

      } else {

        factor[j] = 0.0;

      }

    }

    d_dom_arit[i] = weighted_mean_arit(x, factor);
    d_dom_sqrt[i] = weighted_mean_sqrt(x, factor);
    d_dom_geom[i] = weighted_mean_geom(x, factor);
    d_dom_harm[i] = weighted_mean_harm(x, factor);

    h_dom_arit[i] = weighted_mean_arit(h_0, factor);
    h_dom_sqrt[i] = weighted_mean_sqrt(h_0, factor);
    h_dom_geom[i] = weighted_mean_geom(h_0, factor);
    h_dom_harm[i] = weighted_mean_harm(h_0, factor);

  }

  return DataFrame::create(Named("BAF") = baf_seq,
                           Named("d.0") = d_dom_arit,
                           Named("dg.0") = d_dom_sqrt,
                           Named("dgeom.0") = d_dom_geom,
                           Named("dharm.0") = d_dom_harm,
                           Named("h.0") = h_dom_arit,
                           Named("hg.0") = h_dom_sqrt,
                           Named("hgeom.0") = h_dom_geom,
                           Named("hharm.0") = h_dom_harm);

}
