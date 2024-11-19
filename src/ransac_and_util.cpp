# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]


//-------------------------------------------------------------- "fit_circle_cpp_modified"

// [[Rcpp::export]]
arma::mat fit_circle_cpp_modified(arma::mat points) {
  
  double x1 = points(0, 0);
  double y1 = points(0, 1);
  double x2 = points(1, 0);
  double y2 = points(1, 1);
  double x3 = points(2, 0);
  double y3 = points(2, 1);

  double D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
  double Ux = ((x1 * x1 + y1 * y1) * (y2 - y3) + (x2 * x2 + y2 * y2) * (y3 - y1) + (x3 * x3 + y3 * y3) * (y1 - y2)) / D;
  double Uy = ((x1 * x1 + y1 * y1) * (x3 - x2) + (x2 * x2 + y2 * y2) * (x1 - x3) + (x3 * x3 + y3 * y3) * (x2 - x1)) / D;

  double radius = sqrt((x1 - Ux) * (x1 - Ux) + (y1 - Uy) * (y1 - Uy));
  
  arma::mat center_radius(1,3);
  center_radius.fill(arma::datum::nan);
  
  center_radius(0,0) = Ux;
  center_radius(0,1) = Uy;
  center_radius(0,2) = radius;
  
  return center_radius;
}

//-------------------------------------------------------------- "is_one_row_all_na"

// [[Rcpp::export]]
bool is_one_row_all_na(arma::mat df) {
  return !df.is_finite();         // negate so that in case of TRUE all items of the row will be NA
}


//--------------------------------------------------------------  "sample indices"

// [[Rcpp::export]]
arma::uvec sample_indices(int total_number_obs, int return_obs) {
  return arma::randperm<arma::uvec>(total_number_obs, return_obs);
}

//-------------------------------------------------------------- "internal_ransac"

// [[Rcpp::export]]
arma::field<arma::mat> internal_ransac(arma::mat dat,
                                       double dist = 0.05) {

  unsigned int NROW = dat.n_rows;
  
  arma::uvec idx = sample_indices(NROW, 3);
  arma::uvec idx_cols = {0,1};
  arma::mat inliers = dat(idx, idx_cols);
  
  // First circle fit
  arma::mat fit = fit_circle_cpp_modified(inliers);
  double center_1st = fit(0,0);
  double center_2nd = fit(0,1);
  double radius = fit(0,2);
  
  // Calculate distances for all points and create index
  arma::rowvec radio(NROW);
  radio.fill(arma::datum::nan);
  
  arma::rowvec keep(NROW);
  keep.fill(arma::datum::nan);
  
  int inlier_count = 0;
  
  for(unsigned int i = 0; i < NROW; i++) {

    double radio_i = std::sqrt(std::pow(dat(i,0) - center_1st, 2.0) + std::pow(dat(i,1) - center_2nd, 2.0));
    radio(i) = radio_i;

    // Find inliers
    double abs_val = std::abs(radius - radio_i);
    if (abs_val < dist) {
      keep(i) = i;
      inlier_count++;
    }
  }
  
  arma::mat fit_inlier(1,4);
  fit_inlier.fill(arma::datum::nan);
  fit_inlier(0,0) = center_1st;
  fit_inlier(0,1) = center_2nd;
  fit_inlier(0,2) = radius;
  fit_inlier(0,3) = inlier_count;
  
  // we have to initiallize the matrix outside of the if-condition, otherwise I received :
  // "Error: Mat::operator(): index out of bounds"
  arma::mat df_out(1,6);

  // Check if enough inliers
  if (inlier_count < 2) {
    df_out.fill(arma::datum::nan);
  } else {
    // return data subset
    arma::uvec keep_finite = arma::find_finite(keep);
    
    if (keep_finite.is_empty()) {
      df_out.fill(arma::datum::nan);
    } else {

      df_out.set_size(keep_finite.n_elem, dat.n_cols);
      df_out = dat.rows(keep_finite); 
    }
  }
  
  arma::field<arma::mat> lst_mt(2,1);
  lst_mt(0,0) = df_out;
  lst_mt(1,0) = fit_inlier;

  return lst_mt;
}
  
  
//--------------------------------------------------------------  "RANSAC_cpp"

// [[Rcpp::export]]
arma::mat RANSAC_cpp(arma::mat data,
                     double dist = 0.05) {

  //--------------------------------------------------------------- 1st.
  arma::field<arma::mat> lst_1st = internal_ransac(data, dist);
  arma::mat df_1st = lst_1st(0,0);

  if (is_one_row_all_na(df_1st)) {
    df_1st(0,3) = data.n_rows;
    return df_1st;
  }
  //--------------------------------------------------------------- 2nd.
  arma::field<arma::mat> lst_2nd = internal_ransac(df_1st, dist);
  arma::mat df_2nd = lst_2nd(0,0);

  if (is_one_row_all_na(df_2nd)) {
    df_2nd(0,3) = df_1st.n_rows;
    return df_2nd;
  }

  arma::mat fit = lst_2nd(1,0);
  double center_1st = fit(0,0);
  double center_2nd = fit(0,1);
  double radius = fit(0,2);
  int inlier_count = fit(0,3);
  
  unsigned int NROW = df_2nd.n_rows;
  
  // Calculate final metrics
  double mae = 0.0;
  arma::rowvec distances(NROW);
  distances.fill(arma::datum::nan);
  
  for(unsigned int i = 0; i < NROW; i++) {
    double final_radio = std::sqrt(std::pow(df_2nd(i,0) - center_1st, 2) + std::pow(df_2nd(i,1) - center_2nd, 2));
    mae += (radius - final_radio);
    distances(i) = final_radio;
  }
  
  mae = std::abs(mae) / NROW;
  
  // Calculate cv
  double mean_dist = 0.0;
  double var = 0.0;
  unsigned int N_ELEM = distances.n_elem;
  
  for(unsigned int i = 0; i < N_ELEM; i++) {
    mean_dist += distances(i);
  }
  mean_dist /= N_ELEM;
  
  for(unsigned int i = 0; i < N_ELEM; i++) {
    var += std::pow(distances(i) - mean_dist, 2.0);
  }
  var /= (N_ELEM - 1);
  
  double sd = std::sqrt(var);
  double cv = sd / radius;
  
  // 'x', 'y', 'radio', 'n', 'mae', 'cv'
  arma::mat res_final(1,6);
  res_final.fill(arma::datum::nan);
  
  res_final(0,0) = center_1st;
  res_final(0,1) = center_2nd;
  res_final(0,2) = radius;
  res_final(0,3) = inlier_count;
  res_final(0,4) = mae;
  res_final(0,5) = cv;
  
  return res_final;
}


// function that performs the "RANSAC_cpp" N-times

// [[Rcpp::export]]
arma::mat iterations_RANSAC(arma::mat data, 
                            int n_iterations = 600) {

  // Initialize output matrix
  arma::mat output(n_iterations, 6);                           // 6 columns for: 'x', 'y', 'radio', 'n', 'mae', 'cv'

  for (unsigned int i = 0; i < n_iterations; i++) {

    arma::mat result = RANSAC_cpp(data);
    
    // Fill output matrix
    output(i, 0) = result(0,0);
    output(i, 1) = result(0,1);
    output(i, 2) = result(0,2);
    output(i, 3) = result(0,3);
    output(i, 4) = result(0,4);
    output(i, 5) = result(0,5);
  }
  
  return output;
}


