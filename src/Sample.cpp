#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <rgen.h>
// [[Rcpp::depends(rgen)]]

#include "Distributions.h"

using namespace Rcpp;

void BackwardSample(arma::cube & theta, arma::mat & m, arma::mat & a,
                    arma::cube & C, arma::mat & G, arma::cube & R_inv,
                    const int & n_samples, int & start_slice,
                    const bool & verbose) {
  const int T = theta.n_cols - 1;
  const int p = theta.n_rows;

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
      h_t = m.col(t) + C.slice(t) * G.t() * R_inv.slice(t + 1) *
            (theta.slice(s).col(t + 1) - a.col(t + 1));

      H_t = C.slice(t) - C.slice(t) * G.t() *
            R_inv.slice(t + 1) * G * C.slice(t);

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

void SampleG(arma::mat & G, arma::cube & W_inv, arma::mat & theta,
             arma::mat & Sigma_g_inv, arma::mat & mu_g, const int & p,
             const int & T, const bool Discount = false,
             const double lambda = 0.0) {

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
  V_g = arma::inv_sympd(V_g);
  arma::colvec a_g = kron1.t() * W_tilde_inv *
                     arma::resize(theta.cols(1, T), T * p, 1) +
                     Sigma_g_inv * mu_g;

  arma::mat g = mvnorm(V_g * a_g, V_g);
  G = reshape(g, p, p);
  return;
};

void SampleAR(arma::mat & G, arma::cube & W_inv, arma::mat & theta,
              arma::mat & Sigma_G_inv, arma::mat & mu_G, const int & T,
              const int Discount = false, const double lambda = 0.0) {
  const int p = G.n_rows;
  arma::mat tmp = arma::zeros(p, p);
  arma::mat sum = tmp;
  int W_inv_idx = 0;
  const bool dynamic_W = (W_inv.n_slices > 1);

  for (int t = 1; t <= T; ++t) {
    if (dynamic_W) {
      ++W_inv_idx;
    }
    tmp.diag() = theta.col(t);
    sum += tmp * W_inv.slice(W_inv_idx) * tmp;
  }

  if (Discount) {
    sum *= ((1 + lambda) / lambda);
  }

  arma::mat Sigma_G_new = sum + Sigma_G_inv;
  tmp = arma::solve(Sigma_G_new, sum + Sigma_G_inv * mu_G);

  G.diag() = mvnorm(tmp.diag(), Sigma_G_new);
  //G = W_inv.slice(T);
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

void SampleW_inv (arma::mat & theta, arma::mat & G, arma::mat & W,
                  arma::mat & C_W, const int & df_W, const int & T) {
  arma::mat C_new = theta.cols(1, T) -
                    G * theta.cols(0, T - 1);
  C_new = C_new * C_new.t() + df_W * C_W;
  C_new = arma::inv_sympd(C_new);
  int df_new = df_W + T;
  W = rgen::rwishart(df_new, C_new);
  return;
}
