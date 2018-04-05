#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <rgen.h>
// [[Rcpp::depends(rgen)]]

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

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
  for (int t = 1; t <= T; ++t) {
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
  for (int t = 1; t <= T; ++t) {
    P = G * C.slice(t - 1) * G.t();
    arma::colvec x = theta.slice(i).col(t) - G * theta.slice(i).col(t - 1);
    tmp = (x.t() * solve(P, x));
    total += tmp(0);
  }
  const double beta_new = beta_lambda + total / 2;
  lambda[i + 1] = rigamma(alpha_new, beta_new);
  return;
}

void SampleG(arma::mat & G, arma::mat & W, arma::cube theta,
             arma::mat & Sigma_g_inv, arma::colvec mu_g,
             int & i, const int & p, const int & T, const int S) {

  // FIX ME: Convert to RcppEigen and use selfadjointView
  // https://stackoverflow.com/questions/46700560/converting-an-armadillo-matrix-to-an-eigen-matrixd-and-vice-versa
  // This works if W is fixed
  arma::mat W_inv = arma::inv_sympd(W);
  arma::mat tmp  = theta.slice(i).cols(0, T - 1) *
                   theta.slice(i).cols(0, T - 1).t();
  arma::mat V_g = arma::kron(tmp, W_inv) + Sigma_g_inv;
  Rcout << "trying to invert V_g" << std::endl;
  V_g = arma::inv_sympd(V_g);
  Rcout << "Finished V_g" << std::endl;
  arma::colvec a_g = arma::kron(theta.slice(i).cols(0, T - 1), W_inv)
                     * arma::resize(theta.slice(i).cols(1, T), T * p, 1)
                     + Sigma_g_inv * mu_g;
  Rcout << "Finished a_g" << std::endl;

  // Less efficient approach
  //arma::mat kron1 = kron(theta.slice(i).cols(0, T - 1).t(), arma::eye(p, p));
  //arma::mat W_tilde_inv = arma::kron(arma::eye(T, T), W_inv);
  //Rcout << "Creating V_g" << std::endl;
  //arma::mat V_g = kron1.t() * W_tilde_inv * kron1 + Sigma_g_inv;

  arma::mat g = mvnorm(V_g * a_g, V_g);
  G = reshape(g, p, p);
  return;
};

void SampleV_inv (arma::mat & Y, arma::mat & F, arma::cube & theta,
                  arma::cube & V, arma::mat & C_V, const int & df_V,
                  int & i, const int & T) {
  arma::mat C_new = Y.cols(1, T) - F * theta.slice(i).cols(1, T);
  C_new = C_new * C_new.t() + df_V * C_V;
  C_new = arma::inv_sympd(C_new);
  int df_new = df_V + T;
  V.slice(i) = rgen::rwishart(df_new, C_new);
  return;
}

void SampleW_inv (arma::cube & theta, arma::mat & G,
                  arma::cube & W, arma::mat & C_W, const int & df_W,
                  int & i, const int & T) {
  arma::mat C_new = theta.slice(i).cols(1, T) -
                    G * theta.slice(i).cols(0, T - 1);
  C_new = C_new * C_new.t() + df_W * C_W;
  C_new = arma::inv_sympd(C_new);
  int df_new = df_W + T;
  W.slice(i) = rgen::rwishart(df_new, C_new);
  return;
}
