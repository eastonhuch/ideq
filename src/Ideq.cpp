// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "Kalman.h"
#include "Sample.h"
using namespace Rcpp;

//' Performs FFBS
//'
//' @param Y S by T matrix containing response variable
//' @param F_ S by p matrix defining \eqn{Y_t = F \theta_t + V}
//' @param V S by S variance-covariance matrix of \eqn{V}
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
List Ideq(arma::mat Y, arma::mat F_, arma::mat V,
               arma::mat G, arma::mat W,
               arma::colvec m_0, arma::mat C_0,
               const int n_samples, const bool verbose = false) {
  // This is for known m_0 and C_0
  // figure out dimensions of matrices and check conformability
  const int T = Y.n_cols;
  const int S = Y.n_rows;
  const int p = F_.n_cols;
  CheckDims(Y, F_, V, G, W, m_0, C_0, T, S, p);
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
  arma::cube theta(p, T + 1, n_samples);
  const double alpha_sigma = 0.1;
  const double beta_sigma = 0.1;
  arma::colvec sigma(n_samples);

  for (int i = 0; i < n_samples; i++) {
    if (verbose) {
      Rcout << "Filtering sample number " << i << std::endl;
    }
    checkUserInterrupt();
    Kalman(Y, F_, V, G, W, m, C, T, S, a, R);
    BackwardSample(theta, m, a, C, G, R, T, 1, i, verbose, p);
    SampleSigma(alpha_sigma, beta_sigma, S, T, i, Y, F_, a, sigma);
    // Sample tau^2;
  }

  List out(3);
  out.names() = CharacterVector::create("theta", "sigma2", "tau");
  out[0] = theta;
  out[1] = sigma;
  //out[2] = tau
  return out;
}

// The below R code is for testing
// Simply reload (Ctrl + Shift + L) and create documentation (Ctrl + Shift + D)
/*** R
load('/home/easton/Documents/School/Research/data/test_data.Rdata')
require(fields)
t <- 10; ndraws <- 100
# quilt.plot(latlon[, 1], latlon[, 2], anoms[, 1], ny=20)
small_idx <- latlon[, 1] < 170 & latlon[, 2] > 5
latlon_small <- latlon[small_idx, ]
anoms_small <- anoms[small_idx, 1:t]
quilt.plot(latlon_small[, 1], latlon_small[, 2], anoms_small[, 1], nx = 10, ny = 10)

n <- nrow(anoms_small)
Ft <- Vt <- Gt <- C0 <- diag(n)
Wt <- exp(-as.matrix(dist(latlon_small)))
m0 <- anoms_small[, 1]
dat <- Ideq(anoms_small, Ft, Vt, Gt, Wt, m0, C0, ndraws)

# Plot results compared to raw data
t <- 3 # Time you want to plot
par(mfrow = c(1, 2), mai = c(.4, .5, .2, .2), oma = c(0, 0, 0, .6))
my_breaks <- seq(-1, 1, .1); my_levels <- length(my_breaks) - 1
quilt.plot(latlon_small[, 1], latlon_small[, 2], anoms_small[, t], nx = 10, ny = 10,
           breaks = my_breaks, nlevel = my_levels, add.legend = FALSE)
quilt.plot(latlon_small[, 1], latlon_small[, 2],
           apply(dat[["theta"]][, t + 1 ,], 1, mean), # mean(thetas)
           breaks = my_breaks, nlevel = my_levels,
           nx = 10, ny = 10, ylab = "", yaxt = "n")
# Plot samples of variance parameters
dev.off()
plot(density(dat[["sigma2"]]), xlab = "Sigma2", main = "Sigma2 KDE")
plot(density(dat[["tau"]]), xlab = "Tau", main = "Tau KDE")


? Ideq
*/
