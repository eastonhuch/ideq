// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "Kalman.h"
#include "Sample.h"
#include "Distributions.h"
using namespace Rcpp;

//' Performs FFBS
//'
//' @param Y S by T matrix containing response variable
//' @param F_ S by p matrix defining \eqn{Y_t = F \theta_t + V}
//' @param G p by p matrix defining \eqn{\theta_t = G \theta_{t-1} + W}
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
List Ideq(arma::mat Y, arma::mat F_, arma::mat G,
               arma::colvec m_0, arma::mat C_0,
               const int n_samples, const bool verbose = false) {
  // figure out dimensions of matrices and check conformability
  const int T = Y.n_cols;
  const int S = Y.n_rows;
  const int p = F_.n_cols;
  CheckDims(Y, F_, G, m_0, C_0, T, S, p);
  if (verbose) {
    Rcout << "Dimensions correct, beginning FFBS" << std::endl;
  }

  // create objects for FFBS
  Y.insert_cols(0, 1); // Y is now true-indexed
  arma::mat a(p, T + 1), m(p, T + 1);
  arma::cube R(p, p, T + 1), C(p, p, T + 1);
  m.col(0) = m_0;
  C.slice(0) = C_0;
  arma::cube theta(p, T + 1, n_samples);
  const double alpha_sigma2   = 2.0025,
               alpha_lambda   = 2.25,
               beta_sigma2    = 0.010025,
               beta_lambda    = 0.0625; //FIX ME: these are hard-coded
  arma::colvec sigma2(n_samples + 1), lambda(n_samples + 1);
  // Initial values for sigma2/lambda
  sigma2[0] = rigamma(alpha_sigma2, beta_sigma2);
  lambda[0] = rigamma(alpha_lambda, beta_lambda);

  // FFBS
  for (int i = 0; i < n_samples; ++i) {
    if (verbose) {
      Rcout << "Filtering sample number " << i + 1 << std::endl;
    }
    checkUserInterrupt();
    KalmanDiscounted(Y, F_, G, m, C, a, R, T, S, p, sigma2(i), lambda(i));
    BackwardSample(theta, m, a, C, G, R, T, 1, i, verbose, p);
    SampleSigma2(alpha_sigma2, beta_sigma2, S, T, i, Y, F_, theta, sigma2);
    SampleLambda(alpha_lambda, beta_lambda, p, T, i, G, C , theta, lambda);
  }

  return List::create(_["theta"]  = theta,
                      _["sigma2"] = sigma2,
                      _["lambda"] = lambda);
}

// The below R code is for testing
// Simply reload (Ctrl + Shift + L) and create documentation (Ctrl + Shift + D)
/*** R
# load ocean temperature anomaly data
load('../data/test_data.Rdata')
require(fields)
ts <- 20; ndraws <- 800

# Choose alpha/beta with Method of Moments Estimators
get_prior <- function(m, v) {
  a <- 2 + m^2 / v
  b <- (a - 1) * m
  c(a = a, b = b)
}
get_prior(0.01, 0.2^2) # For sigma2
get_prior(0.05, 0.1^2) # For lambda
# For now these values are hard-coded

# Take small sample of data for debugging
small_idx <- latlon[, 1] < 170 & latlon[, 2] > 5
latlon_small <- latlon[small_idx, ]
anoms_small <- anoms[small_idx, 1:ts]
quilt.plot(latlon_small[, 1], latlon_small[, 2], anoms_small[, 1], nx = 10, ny = 10)

# Create vectors/matrices and fit model
n <- nrow(anoms_small)
Ft <- Gt <- diag(n)
C0 <- exp(-1.5 * as.matrix(dist(latlon_small)))
m0 <- anoms_small[, 1]
dat_full <- Ideq(anoms_small, Ft, Gt, m0, C0, ndraws, verbose = TRUE)
save(dat_full, file = "../data/dat_sample6.RData")
#load("../data/dat_sample2.RData")

# Assess convergence
dev.off()
plot(dat_full[["theta"]][1 ,1 ,], type = "l") # s = 1, t = 1
plot(dat_full[["theta"]][1 ,5 ,], type = "l") # s = 1, t = 5
plot(dat_full[["sigma2"]], type = "l")
plot(dat_full[["lambda"]], type = "l", ylim = c(0, max(dat_full[["lambda"]])))

# lets say the burn-in was 100
burnin <- 300
dat <- list("theta"  = dat_full[["theta"]][, , burnin:ndraws],
            "sigma2" = dat_full[["sigma2"]][burnin:ndraws],
            "lambda" = dat_full[["lambda"]][burnin:ndraws])

# Plot results compared to raw data
par(mfrow = c(1, 2), mai = c(.4, .5, .2, .2), oma = c(0, 0, 0, .6))
my_breaks <- seq(-1.5, 1.5, .1); my_levels <- length(my_breaks) - 1
plot_t <- function(t) {
  quilt.plot(latlon_small[, 1], latlon_small[, 2], anoms_small[, t], nx = 10, ny = 10,
             breaks = my_breaks, nlevel = my_levels, add.legend = FALSE)
  quilt.plot(latlon_small[, 1], latlon_small[, 2],
             apply(dat[["theta"]][, t + 1 ,], 1, mean), # mean(thetas)
             breaks = my_breaks, nlevel = my_levels,
             nx = 10, ny = 10, ylab = "", yaxt = "n", add.legend = FALSE)
}
for (i in 1:ts) {
  plot_t(i)
  Sys.sleep(1)
}

# Plot samples of variance parameters
plot(density(dat[["sigma2"]]), xlab = "Sigma2", main = "Sigma2 KDE")
plot(density(dat[["lambda"]]), xlab = "lambda", main = "lambda KDE")

? Ideq
*/
