#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "Kalman.h"
#include "Sample.h"
#include "Distributions.h"

using namespace Rcpp;

//' Fits a dynamic spatio-temporal model (DSTM) with a discount factor
//'
//' @param Y S by T matrix containing response variable at S spatial locations and T time points
//' @param n_samples integer; number of posterior samples to take
//' @param p integer; dimension of G in the state equation \eqn{\theta_{t+1} = G \theta_{t}}
//' @param verbose boolean; controls verbosity
//' @param sample_sigma2 whether boolean; to sample \eqn{\sigma^2}
//'
//' @keyword IDE, Kalman, Filter
//' @export
//' @examples
//' # Duhh...nothing yet
//' @importFrom Rcpp sourceCpp evalCpp
//' @useDynLib ideq
// [[Rcpp::export]]
List dstm_discount(arma::mat Y, arma::mat F, arma::mat G_0, arma::mat Sigma_G_inv,
                   arma::colvec m_0, arma::mat C_0, NumericVector params,
                   CharacterVector proc_model, const int n_samples,
                   const bool verbose) {
  bool AR = false, FULL = false;
  if (proc_model(0) == "AR") {
    AR = true;
  } else if (proc_model(0) == "Full") {
    FULL = true;
  }

  bool sample_sigma2 = true;
  if (params[4] == NA_REAL) {
    sample_sigma2 = false;
  }

  const int p = G_0.n_rows;
  const int T = Y.n_cols;
  const int S = Y.n_rows;

  // Other objects for sampling
  Y.insert_cols(0, 1); // make Y true-indexed
  arma::cube theta(p, T + 1, n_samples), G;
  arma::mat a(p, T + 1), m(p, T + 1);
  arma::cube R(p, p, T + 1), C(p, p, T + 1), W_inv;
  m.col(0) = m_0;
  C.slice(0) = C_0;

  if (AR || FULL) {
    G.set_size(p, p, n_samples + 1);
    W_inv.set_size(p, p, T + 1);
  } else {
    G.set_size(p, p, 1);
  }
  G.slice(0) = G_0;

  double alpha_lambda = params[0];
  double beta_lambda  = params[1];
  arma::colvec sigma2, lambda(n_samples + 1);
  lambda[0] = rigamma(alpha_lambda, beta_lambda);

  double alpha_sigma2, beta_sigma2, sigma2_i;
  if (sample_sigma2) {
    alpha_sigma2 = params[2];
    beta_sigma2  = params[3];
    sigma2.set_size(n_samples + 1);
    sigma2[0] = rigamma(alpha_sigma2, beta_sigma2);
  } else {
    sigma2_i = params[4];
  }

  int G_idx = 0;
  for (int i = 0; i < n_samples; ++i) {
    if (verbose) {
      Rcout << "Filtering sample number " << i + 1 << std::endl;
    }
    checkUserInterrupt();

    // FFBS
    if (sample_sigma2) sigma2_i = sigma2(i);
    KalmanDiscounted(Y, F, G.slice(G_idx), m, C, a, R, sigma2_i, lambda(i));
    BackwardSample(theta, m, a, C, G.slice(G_idx), R, 1, i, verbose);

    // G
    if (AR) {
      UpdateW_inv(W_inv, C, G.slice(G_idx), AR, lambda(i));
      SampleAR(G.slice(G_idx + 1), W_inv, theta.slice(i), Sigma_G_inv, G_0, T);
    } else if (FULL) {
      UpdateW_inv(W_inv, C, G.slice(G_idx), AR, lambda(i));
      // FIX ME: update W_inv
      //SampleG(G.slice(i + 1), W, theta, Sigma_G_inv, G_0, i, p, T);
    }

    // Sigma2
    if (sample_sigma2) {
      SampleSigma2(alpha_sigma2, beta_sigma2, S, T, i, Y, F, theta, sigma2);
    }

    // Lambda (W)
    SampleLambda(alpha_lambda, beta_lambda, p, T, i, G.slice(G_idx), C, theta, lambda);
    if (AR || FULL) {
      ++ G_idx;
    }
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

//' Fits a dynamic spatio-temporal model (DSTM) that samples the matrix G defining the state equation
//'
//' @param Y S by T matrix containing response variable at S spatial locations and T time points
//' @param n_samples integer; number of posterior samples to take
//' @param p integer; dimension of G in the state equation \eqn{\theta_{t+1} = G \theta_{t}}
//' @param verbose boolean; controls verbosity
//' @param sample_sigma2 whether boolean; to sample \eqn{\sigma^2}
//'
//' @keyword IDE, Kalman, Filter
//' @export
//' @examples
//' # Duhh...nothing yet
//' @importFrom Rcpp sourceCpp evalCpp
//' @useDynLib ideq
// [[Rcpp::export]]
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
    BackwardSample(theta, m, a, C, G.slice(G_idx), R, 1, i, verbose);
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

//' Fits a dynamic spatio-temporal model (DSTM) that samples a diagonal matrix G defining the state equation
//'
//' @param Y S by T matrix containing response variable at S spatial locations and T time points
//' @param n_samples integer; number of posterior samples to take
//' @param p integer; dimension of G in the state equation \eqn{\theta_{t+1} = G \theta_{t}}
//' @param verbose boolean; controls verbosity
//' @param sample_sigma2 whether boolean; to sample \eqn{\sigma^2}
//'
//' @keyword IDE, Kalman, Filter
//' @export
//' @examples
//' # Duhh...nothing yet
//' @importFrom Rcpp sourceCpp evalCpp
//' @useDynLib ideq
// [[Rcpp::export]]
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
    BackwardSample(theta, m, a, C, G.slice(G_idx), R, 1, i, verbose);
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

//' Fits a integrodifference equation model (IDE)
//'
//' @keyword IDE, Kalman, Filter
//' @export
//' @examples
//' # Duhh...nothing yet
//' @importFrom Rcpp sourceCpp evalCpp
//' @useDynLib ideq
// [[Rcpp::export]]
List dstm_IDE() {
  Rcout << "The answer is 42" << std::endl;
  List results;
  results["answer"] = 42;
  return results;
}
