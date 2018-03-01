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
               const int n_samples) {
  // This is for known m_0 and C_0
  // figure out dimensions of matrices and check conformability
  const int T = Y.n_cols;
  const int j = Y.n_rows;
  const int p = F_.n_cols;
  CheckDims(Y, F_, V, G, W, m_0, C_0, T, j, p);
  Rcout << "Dimensions correct, beginning filtering" << std::endl;

  // Kalman Filter
  arma::mat a(p, T);
  arma::cube R(p, p, T);
  arma::mat m(p, T + 1);
  arma::cube C(p, p, T + 1);
  m.col(0) = m_0;
  C.slice(0) = C_0;
  Kalman(Y, F_, V, G, W, m, C, T, j, a, R);
  Rcout << "Finished filtering, beginning sampling" << std::endl;

  // Begin sampling
  arma::colvec h_t(p);
  arma::mat H_t(p, p);
  arma::cube theta(p, T + 1, n_samples);

  for (int s = 0; s < n_samples; s++) {
    theta.slice(s).col(T) = mvnorm(m.col(T), C.slice(T)); // draw theta_T
    // Draw values for theta_{T-1} down to theta_0
    for (int t = T-1; t >= 0; t--) {
      Rcout << "Drawing a value for t = " << t << std::endl;
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
n <- 1600
load('/home/easton/Documents/School/Research/data/test_data.Rdata')
anoms <- anoms[1:n, ]
latlon <- latlon[1:n, ]
Ft <- diag(n)
Gt <- diag(n)
Vt <- .5*diag(n)
D <- as.matrix(dist(latlon))
Wt <- .5*exp(-D^2/10)
m0 <- anoms[,1]
C0 <- diag(n)

quilt.plot(latlon[,1],latlon[,2],anoms[,1],ny=20) # Not working
# This is my work around
quilt.plot_e <- function(x, y, z, title) {
  combined <- cbind(x, y, z)
  x_unique <- sort(unique(latlon[, 1]))
  y_unique <- sort(unique(latlon[, 2]))
  fill_z <- function(vec) {
    a <- combined[combined[, 1] == vec[1] & combined[, 2] == vec[2], 3]
    ifelse(length(a) == 0, NA, a)
  }
  xy <- expand.grid(x_unique, y_unique)
  z <- apply(xy, 1, fill_z)
  z <- matrix(z, nrow = length(x_unique), ncol = length(y_unique))
  image(x_unique, y_unique, z, xlab = "x", ylab = "y", main = title)
}
quilt.plot_e(latlon[, 1], latlon[, 2], anoms[, 1], "Heat Map")

draws <- 1
#dat <- FFBS(anoms, Ft, Vt, Gt, Wt, m0, C0, draws)
#save(dat, file = "/home/easton/Documents/School/Research/data/dat_sample1.RData")
load(file = '/home/easton/Documents/School/Research/data/dat_sample.RData')

plot_FFBS <- function(results, t) {
  par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
  quilt.plot_e(latlon[, 1], latlon[, 2], anoms[, t], "Raw Data")

  results_t0 <- results[, t + 1, ] # t + 1 because index 1 is time 0
  if (draws > 1) results_t0 <- apply(results_t0, 1, mean)
  quilt.plot_e(latlon[, 1], latlon[, 2], results_t0, "FFBS Sample")
  mtext(paste0("Plots for t =", t), outer = TRUE, cex = 1.6)
}

for (i in 1:30) {
  plot_FFBS(dat, i)
  Sys.sleep(1)
}

? FFBS
*/
