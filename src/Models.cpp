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
    G.zeros();
    W_inv.set_size(p, p, T + 1);
  } else {
    G.set_size(p, p, 1);
  }
  G.slice(0) = G_0;

  if (FULL) {
    G_0.reshape(p * p, 1);
  }

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
      CalculateW_inv(W_inv, C, G.slice(G_idx), AR, lambda(i));
      SampleAR(G.slice(G_idx + 1), W_inv, theta.slice(i), Sigma_G_inv, G_0, T);
    } else if (FULL) {
      CalculateW_inv(W_inv, C, G.slice(G_idx), AR, lambda(i));
      SampleG(G.slice(i + 1), W_inv, theta.slice(i), Sigma_G_inv, G_0, p, T);
    }

    // Sigma2
    if (sample_sigma2) {
      SampleSigma2(alpha_sigma2, beta_sigma2, S, T, i, Y, F, theta, sigma2);
    }

    // Is this the right spot to increment G_idx?
    if (AR || FULL) {
      ++ G_idx;
    }

    // Lambda (W)
    SampleLambda(alpha_lambda, beta_lambda, p, T, i,
                 G.slice(G_idx), C, theta, lambda);
  }

  List results;
  results["theta"]  = theta;
  results["lambda"] = lambda;
  if (sample_sigma2) {
    results["sigma2"] = sigma2;
  }
  if (AR || FULL) {
    results["G"] = G;
  }
  results["F"] = F;
  return results;
}

//' Fits a DSTM using a wishart prior for W
//'
//' @keyword Kalman, Filter, FFBS, Wishart
//' @export
//' @examples
//' # Duhh...nothing yet
//' @importFrom Rcpp sourceCpp evalCpp
//' @useDynLib ideq
// [[Rcpp::export]]
List dstm_IW() {
  Rcout << "The answer is 42" << std::endl;
  List results;
  results["answer"] = 42;
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
