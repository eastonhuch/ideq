#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "Kalman.h"
#include "Sample.h"
#include "Distributions.h"

using namespace Rcpp;

List dstm_discount(arma::mat & Y, const int n_samples, const int p,
                   const bool sample_sigma2, const bool verbose) {
  const int T = Y.n_cols;
  const int S = Y.n_rows;

  // Construct F and G
  arma::mat    eig_vecs;
  arma::colvec eig_vals;
  arma::eig_sym(eig_vals, eig_vecs, arma::cov(Y.t()));
  arma::mat F   = eig_vecs.cols(S - p, S - 1);
  arma::mat G_0 = arma::eye(p, p);

  // Other objects for sampling
  Y.insert_cols(0, 1); // make Y true-indexed
  arma::cube theta(p, T + 1, n_samples);
  arma::mat a(p, T + 1), m(p, T + 1);
  arma::cube R(p, p, T + 1), C(p, p, T + 1);
  m.col(0).zeros(); // FIX ME: Create better priors
  C.slice(0).eye();  // FIX ME: Create better priors
  arma::mat W = arma::eye(p, p);

  double alpha_sigma2 = 0, beta_sigma2 = 0,
               alpha_lambda = 0, beta_lambda = 0;
  arma::colvec sigma2, lambda;
  alpha_lambda = 2.25; // FIX ME: How do we want to set priors for sigma?
  beta_lambda  = 0.0625;
  lambda.set_size(n_samples + 1);
  lambda[0] = rigamma(alpha_lambda, beta_lambda);

  double sigma2_i = 0.02; // FIX ME: How do we want to specify this?
  if (sample_sigma2) {
    alpha_sigma2 = 2.0025; // FIX ME: How do we want to set priors for sigma?
    beta_sigma2 = 0.010025;
    sigma2.set_size(n_samples + 1);
    sigma2[0] = rigamma(alpha_sigma2, beta_sigma2);
  }

  for (int i = 0; i < n_samples; ++i) {
    if (verbose) {
      Rcout << "Filtering sample number " << i + 1 << std::endl;
    }
    checkUserInterrupt();

    if (sample_sigma2) sigma2_i = sigma2(i);
    KalmanDiscounted(Y, F, G_0, m, C, a, R, sigma2_i, lambda(i));
    BackwardSample(theta, m, a, C, G_0, R, T, 1, i, verbose, p);

    if (sample_sigma2) {
      SampleSigma2(alpha_sigma2, beta_sigma2, S, T, i, Y, F, theta, sigma2);
    }

    SampleLambda(alpha_lambda, beta_lambda, p, T, i, G_0, C , theta, lambda);
  }

  List results;
  results["theta"]  = theta;
  results["lambda"] = lambda;
  if (sample_sigma2) {
    results["sigma2"] = sigma2;
  }
  results["F"] = F;
  return results;
}

List dstm_sample_G(arma::mat & Y, const int n_samples, const int p,
                   const bool verbose, const bool sample_sigma2) {
  const int T = Y.n_cols;
  const int S = Y.n_rows;

  // Construct F and G
  arma::mat    eig_vecs;
  arma::colvec eig_vals;
  arma::eig_sym(eig_vals, eig_vecs, arma::cov(Y.t()));
  arma::mat F   = eig_vecs.cols(S - p, S - 1);
  arma::cube G(p, p, n_samples + 1);
  G.slice(0).eye();
  arma::colvec mu_g = arma::resize(G.slice(0), p * p, 1);
  arma::mat Sigma_g_inv = arma::eye(p * p, p * p);

  // For sampling sigma2
  double alpha_sigma2 = 0, beta_sigma2 = 0;
  arma::colvec sigma2(1);

  if (sample_sigma2) {
    alpha_sigma2 = 2.0025; // FIX ME: How do we want to set priors for sigma?
    beta_sigma2 = 0.010025;
    sigma2.set_size(n_samples + 1);
    sigma2[0] = rigamma(alpha_sigma2, beta_sigma2);
  } else {
    sigma2[0] = 0.02; // FIX ME: How do we want to specify this?
  }

  // Other objects for sampling
  Y.insert_cols(0, 1); // Y is now true-indexed
  arma::cube theta(p, T + 1, n_samples);
  arma::mat a(p, T + 1), m(p, T + 1);
  arma::cube R(p, p, T + 1), C(p, p, T + 1);
  m.col(0).zeros(); // FIX ME: Create better priors
  C.slice(0).eye();  // FIX ME: Create better priors
  arma::mat W = arma::eye(p, p), V = arma::eye(S, S), V_i; // FIX ME: Make more generic
  if (!sample_sigma2) {
    V_i = sigma2[0] * V;
    V.reset();
  }

  int G_idx = 0; // Avoids complicated control flow
  for (int i = 0; i < n_samples; ++i) {
    if (verbose) {
      Rcout << "Filtering sample number " << i + 1 << std::endl;
    }
    checkUserInterrupt();

    if (sample_sigma2) V_i = sigma2[i] * V;
    Kalman(Y, F, V_i, G.slice(G_idx), W, m, C, a, R);
    BackwardSample(theta, m, a, C, G.slice(G_idx), R, T, 1, i, verbose, p);
    SampleG(G.slice(G_idx + 1), W, theta, Sigma_g_inv, mu_g, i, p, T);
    if (sample_sigma2) {
      SampleSigma2(alpha_sigma2, beta_sigma2, S, T, i, Y, F, theta, sigma2);
    }
    ++G_idx;
  }

  List results;
  results["theta"]  = theta;
  results["G"]      = G;
  results["F"]      = F;
  if (sample_sigma2) {
    results["sigma2"] = sigma2;
  }
  return results;
}

List dstm_AR(arma::mat & Y, const int n_samples, const int p,
             const bool verbose, const bool sample_sigma2) {
  const int T = Y.n_cols;
  const int S = Y.n_rows;

  // Construct F and G
  arma::mat    eig_vecs;
  arma::colvec eig_vals;
  arma::eig_sym(eig_vals, eig_vecs, arma::cov(Y.t()));
  arma::mat F   = eig_vecs.cols(S - p, S - 1);
  arma::cube G = arma::zeros(p, p, n_samples + 1);
  G.slice(0).diag() += 1;
  arma::mat mu_g = G.slice(0);
  arma::mat Sigma_g_inv = arma::eye(p, p);

  // For sampling sigma2
  double alpha_sigma2 = 0, beta_sigma2 = 0;
  arma::colvec sigma2(1);

  if (sample_sigma2) {
    alpha_sigma2 = 2.0025; // FIX ME: How do we want to set priors for sigma?
    beta_sigma2 = 0.010025;
    sigma2.set_size(n_samples + 1);
    sigma2[0] = rigamma(alpha_sigma2, beta_sigma2);
  } else {
    sigma2[0] = 0.02; // FIX ME: How do we want to specify this?
  }

  // Other objects for sampling
  Y.insert_cols(0, 1); // Y is now true-indexed
  arma::cube theta(p, T + 1, n_samples);
  arma::mat a(p, T + 1), m(p, T + 1);
  arma::cube R(p, p, T + 1), C(p, p, T + 1);
  m.col(0).zeros(); // FIX ME: Create better priors
  C.slice(0).eye();  // FIX ME: Create better priors
  arma::mat W = arma::eye(p, p), V = arma::eye(S, S), V_i; // FIX ME: Make more generic
  if (!sample_sigma2) {
    V_i = sigma2[0] * V;
    V.reset();
  }

  int G_idx = 0; // Avoids complicated control flow
  for (int i = 0; i < n_samples; ++i) {
    if (verbose) {
      Rcout << "Filtering sample number " << i + 1 << std::endl;
    }
    checkUserInterrupt();

    if (sample_sigma2) V_i = sigma2[i] * V;
    Kalman(Y, F, V_i, G.slice(G_idx), W, m, C, a, R);
    BackwardSample(theta, m, a, C, G.slice(G_idx), R, T, 1, i, verbose, p);
    SampleAR(G.slice(G_idx + 1), W, theta, Sigma_g_inv, mu_g, i, p, T);
    if (sample_sigma2) {
      SampleSigma2(alpha_sigma2, beta_sigma2, S, T, i, Y, F, theta, sigma2);
    }
    ++G_idx;
  }

  List results;
  results["theta"]  = theta;
  results["G"]      = G;
  results["F"]      = F;
  if (sample_sigma2) {
    results["sigma2"] = sigma2;
  }
  return results;
}

List dstm_IDE() {
  Rcout << "The answer is 42" << std::endl;
  List results;
  results["answer"] = 42;
  return results;
}
