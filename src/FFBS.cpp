// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "Kalman.h"
using namespace Rcpp;

//' Performs FFBS
//'
//' @param Y j by T matrix containing response variable
//' @param F_ j by p matrix defining \eqn{Y_t = F \theta_t + V}
//' @param V j by j variance-covariance matrix of \eqn{V}
//' @param G p by p matrix defining \eqn{\theta_t = G \theta_{t-1} + W}
//' @param W p by p variance-covariance matrix of \eqn{W}
//' @param m_0 p by 1 column vector for a priori mean of \eqn{\theta}
//' @param C_0 p by p matrix of for a priori variance-covariance matrix of \eqn{\theta}
//'
//' @keyword IDE, Kalman, Filter
//' @export
//' @examples
//' # Duhh...nothing yet
//' @useDynLib ideq
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::mat FFBS(arma::mat Y, arma::mat F_, arma::mat V,
               arma::mat G, arma::mat W,
               arma::colvec m_0, arma::mat C_0) {
  // This is for known m_0 and C_0
  Rcout << "chkpt 1" << std::endl;

  // figure out dimensions of matrices and check conformability
  const int T = Y.n_cols;
  const int j = Y.n_rows;
  const int p = F_.n_cols;
  CheckDims(Y, F_, V, G, W, m_0, C_0, T, j, p);
  Rcout << "chkpt 2" << std::endl;

  // Kalman Filter
  arma::mat a(p, T - 1);
  arma::cube R(p, p, T - 1);
  arma::mat m(p, T + 1);
  arma::cube C(p, p, T + 1);
  m.col(0) = m_0;
  C.slice(0) = C_0;
  Rcout << "chkpt 3" << std::endl;
  Kalman(Y, F_, V, G, W, m, C, T, j, a, R);
  Rcout << "chkpt 4" << std::endl;

  // Begin sampling
  arma::colvec h_t(p);
  arma::mat H_t(p, p);
  arma::mat theta(p, T);
  Rcout << "chkpt 5" << std::endl;
  h_t = m.col(T);
  H_t = C.slice(T);
  Rcout << "chkpt 6" << std::endl;

  theta.col(T) = mvnorm(h_t, H_t); // draw theta_T
  Rcout << "chkpt 7" << std::endl;
  // Draw values for theta_{T-1} down to theta_0
  for (int t = T-1; t >= 0; t--) {
    // Mean and variance of theta_t
    h_t = m.col(t) + C.slice(t) * G.t() *
      solve(R.slice(t), (theta.col(t) - a.col(t)));
    Rcout << "chkpt 8" << std::endl;

    H_t = C.slice(t) - C.slice(t) * G.t() *
      solve(R.slice(t), G) * C.slice(t);
    Rcout << "chkpt 9" << std::endl;
    // Draw value for theta_t
    theta.col(t) = mvnorm(h_t, H_t);
    Rcout << "chkpt 10" << std::endl;
  }

  return theta;
}

// The below R code is for testing
// Simply reload (Ctrl + Shift + L) and create documentation (Ctrl + Shift + D)
/*** R
j <- 3
p <- 2
T_ <- 5 # T_ = T
Y <- matrix(c(1.2, 2.4, 2.5, 2.6, 3.2,
             5.6, 7.6, 7.4, 8.1, 8.8,
             0.2, 0.3, 0.3, 0.35, 0.45), nrow = j, ncol = T_, byrow = TRUE)
F_ <- matrix(c(1, 1,
              1, 2.5,
              0.1, 0.2), nrow = j, ncol = p, byrow = TRUE)
V <- diag(c(0.2, 0.5, 0.02))
G <- diag(1, 2)
W <- diag(c(0.1, 0.1))
m_0 <- c(0.3, 0.8)
C_0 <- diag(c(0.1, 0.1))

a <- FFBS(Y, F_, V, G, W, m_0, C_0)
a[[1]] # each column is m_t
a[[2]] # each slice is C_t

? FFBS
*/
