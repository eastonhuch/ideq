#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <rgen.h>
// [[Rcpp::depends(rgen)]]

#include "Distributions.h"
#include "misc_helpers.h"

using namespace Rcpp;

void BackwardSample(arma::mat & theta, const arma::mat & m, const arma::mat & a,
                    const arma::cube & C, const arma::mat & G, const arma::cube & R_inv) {
  const int T = theta.n_cols-1;
  const int p = theta.n_rows;
  arma::colvec h_t(p);
  arma::mat H_t(p, p);

  theta.col(T) = mvnorm(m.col(T), C.slice(T)); // draw theta_T
  // Draw values for theta_{T-1} down to theta_0
  for (int t = T-1; t>=0; --t) {
    // Mean and variance of theta_t
    h_t = m.col(t) + C.slice(t) * G.t() * R_inv.slice(t+1) *
          (theta.col(t+1) - a.col(t+1));

    H_t = C.slice(t) - C.slice(t) * G.t() *
          R_inv.slice(t+1) * G * C.slice(t);
    make_symmetric(H_t);

    // Draw value for theta_t
    theta.col(t) = mvnorm(h_t, H_t);
  }
  return;
}

void SampleSigma2(double & sigma2_new, const double & alpha_sigma2, const double & beta_sigma2,
                  const arma::mat & Y, const arma::mat & F, const arma::mat & theta) {
  const int S = Y.n_rows, T = Y.n_cols-1;
  const double alpha_new = alpha_sigma2 + S*T/2;
  double total = 0;
  for (int t = 1; t <= T; ++t) {
    arma::colvec x = Y.col(t) - F * theta.col(t);
    total += dot(x, x);
  }
  const double beta_new = beta_sigma2 + total / 2;
  sigma2_new = rigamma(alpha_new, beta_new);
  return;
}

void SampleLambda(double & lambda_new, const double & alpha_lambda, const double & beta_lambda,
                  const arma::mat & G, const arma::cube & C, const arma::mat & theta) {
  const int p = G.n_cols, T = theta.n_cols-1;
  const double alpha_new = alpha_lambda + p*T/2;
  arma::mat P(p, p);
  arma::mat tmp(1, 1);
  double total = 0;
  for (int t = 1; t <= T; ++t) {
    P = G * C.slice(t-1) * G.t();
    arma::colvec x = theta.col(t) - G * theta.col(t-1);
    tmp = (x.t() * solve(P, x, arma::solve_opts::equilibrate));
    total += tmp(0);
  }
  const double beta_new = beta_lambda + total / 2;
  lambda_new = rigamma(alpha_new, beta_new);
  return;
}

// NOTE: This update formula comes from Cressie and Wikle p. 457
// Equation (8.54)
void SampleG(arma::mat & G, const arma::cube & W_inv, const arma::mat & theta,
             const arma::mat & Sigma_g_inv, const arma::mat & mu_g,
             const bool Discount = false, const double lambda = 0.0) {
  const int p = theta.n_rows, T = theta.n_cols-1;
  // NOTE: mu_g is arma::reshape(G_0, p * p, 1)
  // Create W_tilde_inv
  arma::mat W_tilde_inv = arma::zeros(T * p, T * p);
  const bool dynamic_W = (W_inv.n_slices > 1);
  int W_idx = 0, start = 0, stop = 0;

  for (int t = 1; t <= T; ++t) {
    if (dynamic_W) {
      ++W_idx;
    }
    start = (t - 1) * p;
    stop  = start + p - 1;
    W_tilde_inv.submat(start, start, stop, stop) = W_inv.slice(W_idx);
  }

  if (Discount) {
    W_tilde_inv *= ((1 + lambda) / lambda);
  }

  // This could probably be optimized
  arma::mat kron1 = kron(theta.cols(0, T - 1).t(), arma::eye(p, p));
  arma::mat V_g = kron1.t() * W_tilde_inv * kron1 + Sigma_g_inv;
  make_symmetric(V_g);
  V_g = arma::inv_sympd(V_g);
  arma::colvec a_g = kron1.t() * W_tilde_inv *
                     arma::reshape(theta.cols(1, T), T * p, 1) +
                     Sigma_g_inv * mu_g;

  arma::mat g = mvnorm(V_g * a_g, V_g);
  G = reshape(g, p, p);
  return;
};

void SampleAR(arma::mat & G, const arma::cube & W_inv, const arma::mat & theta,
              const arma::mat & Sigma_G_inv, const arma::mat & mu_G,
              const bool Discount = false, const double lambda = 0.0) {
  const int T = theta.n_cols-1;
  const int p = G.n_rows;
  arma::mat tmp = arma::zeros(p, p);
  arma::mat sum = tmp;
  arma::colvec sum2 = arma::zeros(p,1);
  int W_inv_idx = 0;
  const bool dynamic_W = (W_inv.n_slices > 1);

  for (int t = 1; t <= T; ++t) {
    if (dynamic_W) {
      ++W_inv_idx;
    }
    tmp.diag() = theta.col(t-1);
    sum  += tmp * W_inv.slice(W_inv_idx) * tmp;
    sum2 += tmp * W_inv.slice(W_inv_idx) * theta.col(t);
  }

  if (Discount) {
    sum  *= ((1 + lambda) / lambda);
    sum2 *= ((1 + lambda) / lambda);
  }

  arma::mat Sigma_G_new = arma::inv_sympd(sum + Sigma_G_inv);
  G.diag() = mvnorm(Sigma_G_new * sum2 + Sigma_G_inv * mu_G, Sigma_G_new);
  return;
};

// Not currently being used
void SampleV(arma::mat & Y, arma::mat & F, arma::mat & theta,
             arma::mat & V, arma::mat & C_V, const int df_V) {
  const int T = theta.n_cols - 1;
  const int p = theta.n_rows;
  arma::mat Y_diffs = Y.cols(1, T) - F * theta.cols(1, T);
  arma::mat C_new = df_V * C_V;
  for (int i=0; i<T; ++i) {
    C_new += Y_diffs.col(i) * Y_diffs.col(i).t();
  }
  int df_new = df_V + T;
  V = rgen::riwishart(df_new, C_new);
  return;
}

void SampleW (arma::mat & W, const arma::mat & theta, const arma::mat & G,
              const arma::mat & C_W, const int df_W) {
  const int T = theta.n_cols - 1;
  const int p = theta.n_rows;
  arma::mat theta_diffs = theta.cols(1, T) -
                          G * theta.cols(0, T - 1);
  arma::mat C_new = df_W * C_W;
  for (int i=0; i<T; ++i) {
    C_new += theta_diffs.col(i) * theta_diffs.col(i).t();
  }
  int df_new = df_W + T;
  W = rgen::riwishart(df_new, C_new);
  return;
}
