#include <RcppCommon.h>
#include <dlib/optimization/optimization.h>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::depends(dlib)]]
//#include <dlib/optimization.h>


using namespace std;
using namespace dlib;


#include <numeric>
#include <vector>
#include <cmath>
#include <complex>

// [[Rcpp::export]]
double norm_euclideanC(arma::vec x) {
  double sq_sum = 0;
  for (int i = 0; i < x.size(); i++) {
    sq_sum += pow(x[i], 2);
  }
  return sqrt(sq_sum);
}

// [[Rcpp::export]]
double dist_euclideanC(arma::vec x, arma::vec y) {
  return norm_euclideanC(x - y);
}

// [[Rcpp::export]]
double dist_euclideanF(arma::vec x, arma::vec y, Function f) {
  NumericVector mod_x = f(x);
  return dist_euclideanC(mod_x, y);
}

// [[Rcpp::export]]
double eta_kernelC(arma::vec t, int lambda) {

  double lambda_num = lambda / 1.0;
  double norm_val = norm_euclideanC(t);
  double y;
  if (lambda % 2 == 0) {
    if (norm_val == 0) {
      y = 0;
    } else {
      y = pow(norm_val, lambda_num) * log(norm_val);
    }
  } else {
    y = pow(norm_val, lambda_num);
  }
  return y;
}

// [[Rcpp::export]]
arma::mat calcE(arma::mat x, int lambda) {
  int nrow = x.n_rows;
  int ncol = x.n_cols;
  arma::mat E(nrow, nrow);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < nrow; j++) {
      E(i, j) = eta_kernelC(x.row(i).t() - x.row(j).t(), lambda);
    }
  }
  return E;
}

// [[Rcpp::export]]
arma::mat etaFunc(arma::vec t, arma::mat tau, int lambda) {
  int ncol = tau.n_cols;
  int nrow = tau.n_rows;
  int vec_length = t.n_elem;

  //Rcpp::Rcout << "ncol: " + std::to_string(ncol) + "\n";
  //Rcpp::Rcout << "nrow: " + std::to_string(nrow) + "\n";
  //Rcpp::Rcout << vec_length;

  arma::mat eta(nrow, 1);

  for (int i = 0; i < nrow; i++) {
    eta(i, 0) = eta_kernelC(t - tau.row(i).t(), lambda);
  }
  return eta;
}

// [[Rcpp::export]]
arma::vec fNew(arma::vec t, arma::mat sol, arma::mat tnew, int I, int d, int lambda) {
  // sol.rows(0, I - 1).print();
  arma::mat temp_mat1 = sol.rows(0, I - 1).t() * etaFunc(t, tnew, lambda);
  // temp_mat1.print();
  arma::vec temp_vec(1);
  temp_vec.ones();
  temp_vec= join_cols(temp_vec, t);
  // temp_vec.print();
  arma::mat temp_mat2 = sol.rows(I, I + d).t() * temp_vec;
  // temp_mat2.print();
  arma::vec out_vec = temp_mat1 + temp_mat2;
  return out_vec;
}

class ProjectionObj {
  public:
    arma::vec y;
    arma::mat sol;
    arma::mat tnew;
    int I;
    int d;
    int lambda;

    double dist_euclidean_opt(dlib::matrix<double, 0, 1> x) {
      std::vector<double> inter(x.begin(), x.end());
      // convert to arma::mat
      arma::vec t = arma::conv_to<arma::vec>::from(inter);
      arma::vec z = fNew(t, sol, tnew, I, d, lambda);
      return dist_euclideanC(y, fNew(t, sol, tnew, I, d, lambda));
    }
};


// [[Rcpp::export]]
arma::vec projectionC(arma::vec x, arma::vec initial_guess, arma::mat sol, arma::mat tnew, int I, int d, int lambda) {
  dlib::matrix<double, 0, 1> init = dlib::mat(initial_guess);
  arma::vec end_points(initial_guess.n_elem);
  ProjectionObj obj;
  obj.y = x;
  obj.sol = sol;
  obj.tnew = tnew;
  obj.I = I;
  obj.d = d;
  obj.lambda = lambda;

  double out = dlib::find_min_using_approximate_derivatives(
    dlib::bfgs_search_strategy(),
    dlib::objective_delta_stop_strategy(1e-5),
    [&obj](dlib::matrix<double, 0, 1> x) {
      return obj.dist_euclidean_opt(x);
    },
    init,
    -1
  );

  for (int i = 0; i < initial_guess.n_elem; i++) {
    end_points[i] = init(i);
  }

  return end_points;
}

// [[Rcpp::export]]
arma::mat calc_tnew(arma::mat x, arma::mat tnew_init, arma::mat sol, int I, int d, int lambda) {
  arma::mat tnew_next(tnew_init.n_rows, tnew_init.n_cols);
  for (int i = 0; i < tnew_next.n_rows; i++) {
    tnew_next.row(i) = projectionC(x.row(i).t(), tnew_init.row(i).t(), sol, tnew_init, I, d, lambda).t();
  }
  return tnew_next;
}

// [[Rcpp::export]]
arma::vec calc_nearest_x(arma::mat df, arma::mat x) {
  arma::vec nearest(df.n_rows);
  for (int i = 0; i < df.n_rows; i++) {
    arma::mat temp_x = x.rows(find(x.col(0) == df(i, 0)));
    arma::vec distances(temp_x.n_rows);
    for (int j = 0; j < temp_x.n_rows; j++) {
      //distances(j) = dist_euclideanC(df.row(i), temp_x.row(j));
      distances(j) = dist_euclideanC(df.row(i).t(), temp_x.row(j).t());
    }
    nearest(i) = distances.index_min() + 1;
  }
  return nearest;
}

