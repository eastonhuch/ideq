#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "Models.h"

using namespace Rcpp;

//' Fits a dynamic spatio-temporal model (DSTM)
//'
//' @param Y S by T matrix containing response variable at S spatial locations and T time points
//' @param F_ S by p matrix defining \eqn{Y_t = F \theta_t + V}
//' @param G_0 p by p matrix defining \eqn{\theta_t = G \theta_{t-1} + W}.
//'        If sample_G is TRUE, then this is used as the starting value and prior mean for G.
//' @param m_0 p by 1 column vector for a priori mean of \eqn{\theta}
//' @param C_0 p by p matrix of for a priori variance-covariance matrix of \eqn{\theta}
//'
//' @keyword IDE, Kalman, Filter
//' @export
//' @examples
//' # Duhh...nothing yet
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
List dstm(arma::mat Y, CharacterVector model = "discount",
          const int n_samples = 1, const int p = 20,
          const bool verbose = false, const bool sample_sigma2 = true) {
  List results;
  if (model(0) == "discount") {
    results = dstm_discount(Y, n_samples, p, sample_sigma2, verbose);
  } else if (model(0) == "sample_G") {
    results = dstm_sample_G(Y, n_samples, p, verbose);
  } else if (model(0) == "IDE") {
    results = dstm_IDE();
  } else if (model(0) == "AR") {
    results = dstm_AR(Y, n_samples, p, verbose);
  } else {
    Rcout << "Model type not recognized" << std::endl;
  }

  return results;
}

// The below R code is for testing
// Simply reload (Ctrl + Shift + L) and create documentation (Ctrl + Shift + D)
/*** R
# load ocean temperature anomaly data
load('../data/test_data.Rdata')
require(fields)
ts <- 30; ndraws <- 200

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

# The lambda values are much larger than expected
dat_full <- dstm(anoms_small, model = "discount", sample_sigma2 = TRUE,
                 n_samples = ndraws, verbose = TRUE)
dat_full <- dstm(anoms, model = "discount", sample_sigma2 = TRUE, p = 8,
                 n_samples = ndraws, verbose = TRUE)
save(dat_full, file = "../data/dat_sample7.RData")
load(file = "../data/dat_sample7.RData")

# This one is working great
dat_full <- dstm(anoms_small, model = "sample_G", sample_sigma2 = TRUE,
                 n_samples = ndraws, verbose = TRUE)

# AR model is doing alright
dat_full <- dstm(anoms_small, model = "AR", sample_sigma2 = TRUE, p = 20,
                 n_samples = ndraws, verbose = TRUE)

# IDE model is still not implemented
dat_full <- dstm(anoms_small, model = "IDE", sample_sigma2 = TRUE,
                 n_samples = ndraws, verbose = TRUE)

# Assess convergence
dev.off()
plot(dat_full[["theta"]][1 ,1 ,], type = "l") # s = 1, t = 1
plot(dat_full[["theta"]][1 ,5 ,], type = "l") # s = 1, t = 5
plot(dat_full[["G"]][1 ,1 ,], type = "l") # s = 1, t = 1
plot(dat_full[["G"]][1 ,10 ,], type = "l") # s = 1, t = 10
plot(dat_full[["sigma2"]], type = "l")
plot(dat_full[["lambda"]], type = "l", ylim = c(0, max(dat_full[["lambda"]])))

# lets say the burn-in was 100
burnin <- 20
dat <- list("theta"  = dat_full[["theta"]][, , burnin:ndraws],
            "sigma2" = dat_full[["sigma2"]][burnin:ndraws],
            "lambda" = dat_full[["lambda"]][burnin:ndraws],
            "F"      = dat_full[["F"]])

# Plot results compared to raw data
par(mfrow = c(1, 2), mai = c(.4, .5, .2, .2), oma = c(0, 0, 0, .6))
my_breaks <- seq(-1.5, 1.5, .1); my_levels <- length(my_breaks) - 1

plot_t <- function(t) {
  quilt.plot(latlon_small[, 1], latlon_small[, 2], anoms_small[, t], nx = 10, ny = 10,
             breaks = my_breaks, nlevel = my_levels, add.legend = FALSE)
  quilt.plot(latlon_small[, 1], latlon_small[, 2],
             dat[["F"]] %*% apply(dat[["theta"]][, t + 1,], 1, mean), # mean(thetas)
             breaks = my_breaks, nlevel = my_levels,
             nx = 10, ny = 10, ylab = "", yaxt = "n", add.legend = FALSE)
}
for (i in 1:ts) {
  plot_t(i)
  Sys.sleep(1)
}

# For full data
plot_t_full <- function(t) {
  quilt.plot(latlon[, 1], latlon[, 2], anoms[, t], nx = 71, ny = 25,
             breaks = my_breaks, nlevel = my_levels, add.legend = FALSE)
  quilt.plot(latlon[, 1], latlon[, 2],
             dat[["F"]] %*% apply(dat[["theta"]][, t + 1,], 1, mean), # mean(thetas)
             breaks = my_breaks, nlevel = my_levels,
             nx = 71, ny = 25, ylab = "", yaxt = "n", add.legend = FALSE)
}
for (i in 1:ts) {
  plot_t_full(i)
  Sys.sleep(1)
}

# Plot samples of variance parameters
plot(density(dat[["sigma2"]]), xlab = "Sigma2", main = "Sigma2 KDE")
plot(density(dat[["lambda"]]), xlab = "lambda", main = "lambda KDE")

? dstm
*/
