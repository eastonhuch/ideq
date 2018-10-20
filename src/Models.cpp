#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "Kalman.h"
#include "Sample.h"
#include "Distributions.h"
#include "IDE_helpers.h"

using namespace Rcpp;

//' Fits a dynamic spatio-temporal model (DSTM) with a discount factor
//'
//' @param Y S by T matrix containing response variable at S spatial locations and T time points
//' @param n_samples integer; number of posterior samples to take
//' @param p integer; dimension of G in the state equation \eqn{\theta_{t+1} = G \theta_{t}}
//' @param verbose boolean; controls verbosity
//' @param sample_sigma2 whether boolean; to sample \eqn{\sigma^2}
//'
//' @export
//' @examples
//' # Duhh...nothing yet
//' @importFrom Rcpp sourceCpp evalCpp
//' @useDynLib ideq
// [[Rcpp::export]]
List dstm_discount(arma::mat Y, arma::mat F, arma::mat G_0, arma::mat Sigma_G_inv,
                   arma::colvec m_0, arma::mat C_0, NumericVector params,
                   CharacterVector proc_model, const int n_samples, const bool verbose) {
  // Create high-level model parameters
  bool AR = proc_model(0) == "AR";
  bool FULL = proc_model(0) == "Full";
  bool sample_sigma2 = params[4] != NA_LOGICAL;
  const int p = G_0.n_rows;
  const int T = Y.n_cols;
  const int S = Y.n_rows;

  // Create matrices and cubes for FFBS
  Y.insert_cols(0, 1); // make Y true-indexed; i.e. index 1 is t_1
  arma::cube theta(p, T + 1, n_samples), G;
  arma::mat a(p, T + 1), m(p, T + 1), tmp;
  arma::cube R_inv(p, p, T + 1), C(p, p, T + 1), C_T(p, p, n_samples+1);
  m.col(0) = m_0;
  C.slice(0) = C_0;

  if (AR) {
    G.set_size(p, p, n_samples+1);
    G.zeros();
    G.slice(0).diag() = mvnorm(G_0.diag(), arma::inv_sympd(Sigma_G_inv));
  } else if (FULL) {
    G.set_size(p, p, n_samples+1);
    G_0.reshape(p*p, 1);
    tmp = mvnorm(G_0, arma::inv_sympd(Sigma_G_inv));
    tmp.reshape(p, p);
    G.slice(0) = tmp;
  } else {
    G.set_size(p, p, 1);
    G.slice(0) = G_0;
  }

  // Create parameters for sampling lambda and sigma2
  double alpha_lambda = params[0];
  double beta_lambda  = params[1];
  arma::vec sigma2, lambda(n_samples + 1);
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

  // Begin MCMC
  int G_idx = 0; // This value is incremented each iteration for AR and Full models
  for (int i=0; i<n_samples; ++i) {
    checkUserInterrupt();

    // FFBS
    if (verbose) {
      Rcout << "Filtering sample number " << i + 1 << std::endl;
    }
    if (sample_sigma2) sigma2_i = sigma2(i);
    KalmanDiscount(m, C, a, R_inv, Y, F, G.slice(G_idx), sigma2_i, lambda(i));
    C_T.slice(i+1) = C.slice(T); // Save for predictions

    if (verbose) {
      Rcout << "Drawing sample number " << i + 1 << std::endl;
    }
    BackwardSample(theta.slice(i), m, a, C, G.slice(G_idx), R_inv);

    // Sigma2
    if (sample_sigma2) {
      SampleSigma2(sigma2(i+1), alpha_sigma2, beta_sigma2, Y, F, theta.slice(i));
    }

    // Lambda (W)
    SampleLambda(lambda(i+1), alpha_lambda, beta_lambda,
                 G.slice(G_idx), C, theta.slice(i));

    // G
    if (AR) {
      SampleAR(G.slice(G_idx+1), R_inv, theta.slice(i),
               Sigma_G_inv, G_0, true, lambda(i+1));
      ++G_idx;
    } else if (FULL) {
      SampleG(G.slice(G_idx+1), R_inv, theta.slice(i),
              Sigma_G_inv, G_0, true, lambda(i+1));
      ++G_idx;
    }
  }

  List results;
  results["theta"]  = theta;
  results["lambda"] = lambda;
  if (sample_sigma2) {
    results["sigma2"] = sigma2;
  } else {
    results["sigma2"] = sigma2_i;
  }
  if (AR || FULL) {
    results["G"] = G;
  }
  results["F"] = F;
  results["C_T"] = C_T;
  return results;
}

//' Fits a DSTM using a wishart prior for W
//'
//' @export
//' @examples
//' # Duhh...nothing yet
//' @importFrom Rcpp sourceCpp evalCpp
//' @useDynLib ideq
// [[Rcpp::export]]
List dstm_IW(arma::mat Y, arma::mat F, arma::mat G_0, arma::mat Sigma_G_inv,
             arma::colvec m_0, arma::mat C_0, arma::mat C_W,
             NumericVector params, CharacterVector proc_model,
             const int n_samples, const bool verbose) {
  // Create high-level model parameters
  bool AR = proc_model(0) == "AR";
  bool FULL = proc_model(0) == "Full";
  const int p = G_0.n_rows;
  const int T = Y.n_cols;
  const int S = Y.n_rows;

  // Create variance parameters
  arma::vec sigma2;
  double alpha_sigma2, beta_sigma2, sigma2_i;
  bool sample_sigma2;
  if (params[3] == NA_LOGICAL) {
    sample_sigma2 = false;
    sigma2_i = params[3];
  }
  else {
    sample_sigma2 = true;
    alpha_sigma2 = params[1];
    beta_sigma2  = params[2];
    sigma2.set_size(n_samples + 1);
    sigma2[0] = rigamma(alpha_sigma2, beta_sigma2);
  }
  const double df_W = params[0];

  // Create matrices and cubes for FFBS
  Y.insert_cols(0, 1); // make Y true-indexed; i.e. index 1 is t_1
  arma::cube theta(p, T + 1, n_samples), G;
  arma::mat a(p, T + 1), m(p, T + 1), tmp;
  arma::cube R_inv(p, p, T + 1), C(p, p, T + 1), W(p, p, n_samples + 1);
  m.col(0) = m_0;
  C.slice(0) = C_0;
  W.slice(0) = df_W * C_W;

  if (AR) {
    G.set_size(p, p, n_samples + 1);
    G.zeros();
    G.slice(0).diag() = mvnorm(G_0.diag(), arma::inv_sympd(Sigma_G_inv));
  } else if (FULL) {
    G.set_size(p, p, n_samples + 1);
    G_0.reshape(p*p, 1);
    tmp = mvnorm(G_0, arma::inv_sympd(Sigma_G_inv));
    tmp.reshape(p, p);
    G.slice(0) = tmp;
  } else {
    G.set_size(p, p, 1);
    G.slice(0) = G_0;
  }

  // Begin MCMC
  int G_idx = 0; // This value is incremented each iteration for AR and Full models
  for (int i = 0; i < n_samples; ++i) {
    checkUserInterrupt();

    // FFBS
    if (verbose) {
      Rcout << "Filtering sample number " << i + 1 << std::endl;
    }
    if (sample_sigma2) sigma2_i = sigma2(i);
    Kalman(m, C, a, R_inv, Y, F, G.slice(G_idx), W.slice(i), sigma2_i);

    if (verbose) {
      Rcout << "Drawing sample number " << i + 1 << std::endl;
    }
    BackwardSample(theta.slice(i), m, a, C, G.slice(G_idx), R_inv);

    // Sigma2
    if (sample_sigma2) {
      SampleSigma2(sigma2(i+1), alpha_sigma2, beta_sigma2, Y, F, theta.slice(i));
    }

    // W
    SampleW(W.slice(i+1), theta.slice(i), G.slice(G_idx), C_W, df_W);

    // G
    if (AR) {
      SampleAR(G.slice(G_idx+1), W.slices(i+1, i+1), theta.slice(i), Sigma_G_inv, G_0);
      ++G_idx;
    } else if (FULL) {
      SampleG(G.slice(G_idx+1), W.slices(i+1, i+1), theta.slice(i), Sigma_G_inv, G_0);
      ++G_idx;
    }

  }

  List results;
  results["theta"]  = theta;
  if (sample_sigma2) {
    results["sigma2"] = sigma2;
  } else {
    results["sigma2"] = sigma2_i;
  }
  if (AR || FULL) {
    results["G"] = G;
  }
  results["F"] = F;
  results["W"] = W;
  return results;
}

//' Fits a integrodifference equation model (IDE)
//'
//' @export
//' @examples @importFrom Rcpp sour
//' # Duhh...nothing yet
//'ceCpp evalCpp
//' @useDynLib ideq
// [[Rcpp::export]]
List dstm_IDE(arma::mat Y, arma::mat locs, arma::colvec m_0, arma::mat C_0,
              NumericVector params, const int n_samples, const bool verbose) {
  // Extract scalar parameters
  const double J  = params["J"];
  const double L  = params["L"];
  const double alpha_sigma2  = params["alpha_sigma2"];
  const double beta_sigma2   = params["beta_sigma2"];
  double       sigma2_i      = params["sigma2"];
  const bool   sample_sigma2 = sigma2_i == NA_LOGICAL;
  const double alpha_lambda  = params["alpha_lambda"];
  const double beta_lambda   = params["alpha_lambda"];
  const int locs_dim = locs.n_cols;
  const int p = 2*J*J + 1;
  const int T = Y.n_cols;
  const int S = Y.n_rows;

  // Create matrices and cubes for FFBS
  Y.insert_cols(0, 1); // make Y true-indexed; i.e. index 1 is t_1
  arma::cube theta(p, T + 1, n_samples), G(p, p, n_samples+1);
  arma::mat a(p, T + 1), m(locs_dim, T + 1);
  arma::cube R_inv(p, p, T + 1), C(locs_dim, locs_dim, T + 1), C_T(p, p, n_samples+1);
  m.col(0) = m_0;
  C.slice(0) = C_0;

  // Create observation matrix (F) and initial process matrix (G)
  arma::mat w_for_B = make_w_for_B(locs, J, L);
  arma::mat F = makeF(locs, w_for_B, J, L);
  arma::mat FtFiFt = arma::solve(F.t() * F, F.t());
  arma::mat B = makeB(m_0, C_0, locs, w_for_B, J, L);
  G.slice(0) = FtFiFt * B;

  // Create variance parameters
  arma::vec sigma2, lambda(n_samples + 1);
  lambda[0] = rigamma(alpha_lambda, beta_lambda);
  if (sample_sigma2) {
    sigma2.set_size(n_samples+1);
    sigma2[0] = rigamma(alpha_sigma2, beta_sigma2);
  }

  // Begin MCMC
  for (int i = 0; i < n_samples; ++i) {
    checkUserInterrupt();

    // FFBS
    if (verbose) {
      Rcout << "Filtering sample number " << i + 1 << std::endl;
    }
    if (sample_sigma2) sigma2_i = sigma2(i);
    KalmanDiscount(m, C, a, R_inv, Y, F, G.slice(i), sigma2_i, lambda(i));

    if (verbose) {
      Rcout << "Drawing sample number " << i + 1 << std::endl;
    }
    BackwardSample(theta.slice(i), m, a, C, G.slice(i), R_inv);
    C_T.slice(i+1) = C.slice(T); // Save for predictions

    // Sigma2
    if (sample_sigma2) {
      SampleSigma2(sigma2(i+1), alpha_sigma2, beta_sigma2, Y, F, theta.slice(i));
    }

    // Lambda (W)
    SampleLambda(lambda(i+1), alpha_lambda, beta_lambda,
                 G.slice(i), C, theta.slice(i));

    // MH step for mu
    // Sample proposal value
    // Create new B and G matrix from value
    // Likelihood looks something like this
    // (theta_t - G * theta_{t-1})' (lambda * C)^{-1} (theta_t - G * theta_{t-1})
    // Where C comes from the filtering recursions

    // MH step for Sigma
    // Similar to portion for mu

    // Make new process matrix
    //arma::mat B = makeB(m_0, C_0, locs, w_for_B, J, L);
    G.slice(i) = FtFiFt * B;
  }

  Rcout << "The answer is 42" << std::endl;
  List results;
  results["answer"] = 42;
  return results;
}
