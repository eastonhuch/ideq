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
arma::cube FFBS(arma::mat Y, arma::mat F_, arma::mat V,
               arma::mat G, arma::mat W,
               arma::colvec m_0, arma::mat C_0,
               const int n_samples, const bool verbose = false) {
  // This is for known m_0 and C_0
  // figure out dimensions of matrices and check conformability
  const int T = Y.n_cols;
  const int j = Y.n_rows;
  const int p = F_.n_cols;
  CheckDims(Y, F_, V, G, W, m_0, C_0, T, j, p);
  if (verbose) {
  Rcout << "Dimensions correct, beginning filtering" << std::endl;
  }

  // Kalman Filter
  arma::mat a(p, T);
  arma::cube R(p, p, T);
  arma::mat m(p, T + 1);
  arma::cube C(p, p, T + 1);
  m.col(0) = m_0;
  C.slice(0) = C_0;
  Kalman(Y, F_, V, G, W, m, C, T, j, a, R);
  if (verbose) {
  Rcout << "Finished filtering, beginning sampling" << std::endl;
  }

  // Begin sampling
  arma::colvec h_t(p);
  arma::mat H_t(p, p);
  arma::cube theta(p, T + 1, n_samples);

  for (int s = 0; s < n_samples; s++) {
    theta.slice(s).col(T) = mvnorm(m.col(T), C.slice(T)); // draw theta_T
    // Draw values for theta_{T-1} down to theta_0
    for (int t = T-1; t >= 0; t--) {
      if (verbose) {
      Rcout << "Drawing a value for t = " << t << std::endl;
      }
      checkUserInterrupt();
      // Mean and variance of theta_t
      h_t = m.col(t) + C.slice(t) * G.t() *
        solve(R.slice(t), (theta.slice(s).col(t + 1) - a.col(t)));

      H_t = C.slice(t) - C.slice(t) * G.t() *
        solve(R.slice(t), G) * C.slice(t);
      // Draw value for theta_t
      theta.slice(s).col(t) = mvnorm(h_t, H_t);
    }
  }
  return theta;
}

// The below R code is for testing
// Simply reload (Ctrl + Shift + L) and create documentation (Ctrl + Shift + D)
/*** R
n <- 1600 # number of locations
p <- 200  # number of basis functions
t <- 10 # number of time points
ndraws <- 1
load('/home/easton/Documents/School/Research/data/test_data.Rdata')
require(fields)
anoms_i <- anoms[1:n, 1:t]
latlon_i <- latlon[1:n, ]
quilt.plot(latlon_i[, 1], latlon_i[, 2], anoms_i[, 1], ny=20)

C <- cov(t(anoms_i))
EOFs <- eigen(C)$vec
Ft <- EOFs[, 1:p]
Vt <- .5 * diag(n)

Gt <- diag(p)
D <- as.matrix(dist(latlon_i))
Wt <- t(Ft) %*% exp(-D^2) %*% Ft

m0 <- matrix(0, nrow = p, ncol = 1)
C0 <- diag(p)

dat <- FFBS(anoms, Ft, Vt, Gt, Wt, m0, C0, ndraws)
#save(dat, file = "/home/easton/Documents/School/Research/data/dat_sample1.RData")
#load(file = '/home/easton/Documents/School/Research/data/dat_sample.RData')

? FFBS
*/
