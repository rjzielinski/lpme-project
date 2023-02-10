#include <Rcpp.h>
using namespace Rcpp;

#include <numeric>
#include <vector>
#include <cmath>
#include <complex>

// [[Rcpp::export]]
double norm_euclideanC(NumericVector x) {
  double sq_sum = 0;
  for (int i = 0; i < x.size(); i++) {
    sq_sum += pow(x[i], 2);
  }
  return sqrt(sq_sum);
}

// [[Rcpp::export]]
double dist_euclideanC(NumericVector x, NumericVector y) {
  return norm_euclideanC(x - y);
}

// [[Rcpp::export]]
double eta_kernelC(NumericVector t, int lambda) {
  double norm_val = norm_euclideanC(t);
  double y;
  if (lambda % 2 == 0) {
    if (norm_val == 0) {
      y = 0;
    } else {
      y = pow(norm_val, lambda) * log(norm_val);
    }
  } else {
    y = pow(norm_val, lambda);
  }
  return y;
}

// [[Rcpp::export]]
NumericMatrix calcE(NumericMatrix x, int lambda) {
  int nrow = x.nrow();
  int ncol = x.ncol();
  NumericMatrix E(nrow, nrow);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < nrow; j++) {
      NumericVector temp_row_i(ncol);
      NumericVector temp_row_j(ncol);
      for (int k = 0; k < ncol; k++) {
        temp_row_i[k] = x(i, k);
        temp_row_j[k] = x(j, k);
      }
      E(i, j) = eta_kernelC(temp_row_i - temp_row_j, lambda);
    }
  }
  return E;
}

// [[Rcpp::export]]
NumericMatrix etaFunc(NumericVector t, NumericMatrix tau, int lambda) {
  int ncol = tau.ncol();
  int nrow = tau.nrow();
  int vec_length = t.size();

  //Rcpp::Rcout << "ncol: " + std::to_string(ncol) + "\n";
  //Rcpp::Rcout << "nrow: " + std::to_string(nrow) + "\n";
  //Rcpp::Rcout << vec_length;

  NumericMatrix eta(nrow, 1);

  for (int i = 0; i < nrow; i++) {
    NumericVector temp_row(ncol);
    for (int j = 0; j < ncol; j++) {
      temp_row[j] = tau(i, j);
    }
    eta(i, 0) = eta_kernelC(t - temp_row, lambda);
  }
  return eta;
}

