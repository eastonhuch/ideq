// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "Distributions.h"
using namespace Rcpp;

void BackwardSample(arma::cube & theta, arma::mat & m, arma::mat & a,
                    arma::cube & C, arma::mat & G,
                    arma::cube & R, const int & T,
                    const int & n_samples, int & start_slice,
                    const bool & verbose, const int & p) {
  arma::colvec h_t(p);
  arma::mat H_t(p, p);

  // NOTE: n_samples is how many samples we want to draw on this function call
  for (int s = start_slice; s < (start_slice + n_samples); ++s) {
    if (verbose) {
      Rcout << "Drawing sample number " << s + 1 << std::endl;
    }
    checkUserInterrupt();
    theta.slice(s).col(T) = mvnorm(m.col(T), C.slice(T)); // draw theta_T
    // Draw values for theta_{T-1} down to theta_0
    for (int t = T-1; t >= 0; --t) {
      // Mean and variance of theta_t
      h_t = m.col(t) + C.slice(t) * G.t() *
        solve(R.slice(t + 1), (theta.slice(s).col(t + 1) - a.col(t + 1)));

      H_t = C.slice(t) - C.slice(t) * G.t() *
        solve(R.slice(t + 1), G) * C.slice(t);
      // Draw value for theta_t
      theta.slice(s).col(t) = mvnorm(h_t, H_t);
    }
  }
  return;
}

void SampleSigma2(const double & alpha_sigma2, const double & beta_sigma2,
                 const int & S, const int & T, int i,
                 arma::mat & Y, arma::mat & F_,
                 arma::cube & theta, arma::colvec & sigma2) {
  const double alpha_new = alpha_sigma2 + S * T / 2;
  double total = 0;
  for (int t = 1; t < T; ++t) {
    arma::colvec x = Y.col(t) - F_ * theta.slice(i).col(t);
    total += dot(x, x);
  }
  const double beta_new = beta_sigma2 + total / 2;
  sigma2[i + 1] = rigamma(alpha_new, beta_new);
  return;
}

void SampleLambda(const double & alpha_lambda, const double & beta_lambda,
                 const int & p, const int & T, int i,
                 arma::mat & G, arma::cube & C,
                 arma::cube & theta, arma::colvec & lambda) {
  const double alpha_new = alpha_lambda + p * T / 2;
  arma::mat P(p, p);
  arma::mat tmp(1, 1);
  double total = 0;
  for (int t = 1; t < T; ++t) {
    P = G * C.slice(t - 1) * G.t();
    arma::colvec x = theta.slice(i).col(t) - G * theta.slice(i).col(t - 1);
    tmp = (x.t() * solve(P, x));
    total += tmp(0);
  }
  const double beta_new = beta_lambda + total / 2;
  lambda[i + 1] = rigamma(alpha_new, beta_new);
  return;
}
