#include <Rcpp.h>
using namespace Rcpp;

#include <numeric>
#include <vector>
#include <cmath>

// [[Rcpp::export]]
double norm_euclidean(NumericVector x) {
  if (x.length() > 1) {
    return sqrt(std::inner_product(x.begin(), x.end(), x.begin(), 0));
  } else {
    return sqrt(pow(x[0], 2));
  }
}

// [[Rcpp::export]]
double dist_euclidean(NumericVector x, NumericVector y) {
  return norm_euclidean(x - y);
}

// [[Rcpp::export]]
double eta_kernel(NumericVector t, int lambda) {
  double norm_val = norm_euclidean(t);
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
