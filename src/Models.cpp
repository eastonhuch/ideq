#include <math.h>
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
  // Extract scalar parameters
  bool AR = proc_model(0) == "AR";
  bool FULL = proc_model(0) == "Full";
  const int p = G_0.n_rows;
  const int T = Y.n_cols;
  const int S = Y.n_rows;
  const double alpha_sigma2 = params["alpha_sigma2"];
  const double beta_sigma2  = params["beta_sigma2"];
  const bool sample_sigma2  = params["sample_sigma2"] > 0;
  double alpha_lambda = params["alpha_lambda"];
  double beta_lambda  = params["beta_lambda"];

  // Create matrices and cubes for FFBS
  Y.insert_cols(0, 1); // make Y true-indexed; i.e. index 1 is t_1
  arma::cube theta(p, T+1, n_samples), G;
  theta.slice(0).zeros();
  arma::mat a(p, T+1), m(p, T+1), tmp;
  arma::cube R_inv(p, p, T+1), C(p, p, T+1), C_T(p, p, n_samples+1);
  m.col(0) = m_0;
  C.slice(0) = C_0;

  if (AR) {
    G.set_size(p, p, n_samples+1);
    G.zeros();
    tmp = G_0.diag();
    G_0.set_size(p, 1);
    G_0 = tmp;
    G.slice(0).diag() = mvnorm(G_0, arma::inv_sympd(Sigma_G_inv));
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

  // Create variance parameters
  arma::vec sigma2;
  double sigma2_i;
  if (sample_sigma2) {
    sigma2.set_size(n_samples+1);
    SampleSigma2(sigma2.at(0), alpha_sigma2, beta_sigma2, Y, F, theta.slice(0));
  }
  else {
    sigma2_i = params["sigma2"];
  }

  // Create parameters for sampling lambda
  arma::vec lambda(n_samples+1);
  lambda.at(0) = rigamma(alpha_lambda, beta_lambda);

  // Begin MCMC
  int G_idx = 0; // This value is incremented each iteration for AR and Full models
  for (int i=0; i<n_samples; ++i) {
    checkUserInterrupt();

    // FFBS
    if (verbose) {
      Rcout << "Filtering sample number " << i+1 << std::endl;
    }
    if (sample_sigma2) sigma2_i = sigma2.at(i);
    KalmanDiscount(m, C, a, R_inv, Y, F, G.slice(G_idx), sigma2_i, lambda.at(i));
    C_T.slice(i+1) = C.slice(T); // Save for predictions

    if (verbose) {
      Rcout << "Drawing sample number " << i+1 << std::endl;
    }
    BackwardSample(theta.slice(i), m, a, C, G.slice(G_idx), R_inv);

    // Sigma2
    if (sample_sigma2) {
      SampleSigma2(sigma2.at(i+1), alpha_sigma2, beta_sigma2, Y, F, theta.slice(i));
    }

    // Lambda (W)
    SampleLambda(lambda.at(i+1), alpha_lambda, beta_lambda,
                 G.slice(G_idx), C, theta.slice(i));

    // G
    if (AR) {
      SampleAR(G.slice(G_idx+1), R_inv, theta.slice(i),
               Sigma_G_inv, G_0, true, lambda.at(i+1));
      ++G_idx;
    } else if (FULL) {
      SampleG(G.slice(G_idx+1), R_inv, theta.slice(i),
              Sigma_G_inv, G_0, true, lambda.at(i+1));
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
  // Extract scalar parameters
  bool AR = proc_model(0) == "AR";
  bool FULL = proc_model(0) == "Full";
  const int p = G_0.n_rows;
  const int T = Y.n_cols;
  const int S = Y.n_rows;
  const double alpha_sigma2 = params["alpha_sigma2"];
  const double beta_sigma2  = params["beta_sigma2"];
  const bool sample_sigma2  = params["sample_sigma2"] > 0;

  // Create matrices and cubes for FFBS
  Y.insert_cols(0, 1); // make Y true-indexed; i.e. index 1 is t_1
  arma::cube theta(p, T+1, n_samples), G;
  theta.slice(0).zeros();
  arma::mat a(p, T+1), m(p, T+1), tmp;
  arma::cube R_inv(p, p, T+1), C(p, p, T+1), W(p, p, n_samples+1);
  m.col(0) = m_0;
  C.slice(0) = C_0;
  const double df_W = params["df_W"];
  W.slice(0) = df_W * C_W;

  if (AR) {
    G.set_size(p, p, n_samples+1);
    G.zeros();
    tmp = G_0.diag();
    G_0.set_size(p, 1);
    G_0 = tmp;
    G.slice(0).diag() = mvnorm(G_0, arma::inv_sympd(Sigma_G_inv));
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

  // Create variance parameters
  arma::vec sigma2;
  double sigma2_i;
  if (sample_sigma2) {
    sigma2.set_size(n_samples+1);
    SampleSigma2(sigma2.at(0), alpha_sigma2, beta_sigma2, Y, F, theta.slice(0));
  }
  else {
    sigma2_i = params["sigma2"];
  }

  // Begin MCMC
  int G_idx = 0; // This value is incremented each iteration for AR and Full models
  for (int i = 0; i < n_samples; ++i) {
    checkUserInterrupt();

    // FFBS
    if (verbose) {
      Rcout << "Filtering sample number " << i+1 << std::endl;
    }
    if (sample_sigma2) sigma2_i = sigma2.at(i);
    Kalman(m, C, a, R_inv, Y, F, G.slice(G_idx), W.slice(i), sigma2_i);

    if (verbose) {
      Rcout << "Drawing sample number " << i+1 << std::endl;
    }
    BackwardSample(theta.slice(i), m, a, C, G.slice(G_idx), R_inv);

    // Sigma2
    if (sample_sigma2) {
      SampleSigma2(sigma2.at(i+1), alpha_sigma2, beta_sigma2, Y, F, theta.slice(i));
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
  if (sample_sigma2) {    results["sigma2"] = sigma2;
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
              arma::colvec m_kernel, arma::mat C_kernel,
              NumericVector params, const int n_samples, const bool verbose) {
  // Extract scalar parameters
  const double J  = params["J"];
  const double L  = params["L"];
  const double alpha_sigma2  = params["alpha_sigma2"];
  const double beta_sigma2   = params["beta_sigma2"];
  double       sigma2_i      = params["sigma2"];
  const bool   sample_sigma2 = sigma2_i == NA;
  const double alpha_lambda  = params["alpha_lambda"];
  const double beta_lambda   = params["beta_lambda"];
  const double proposal_factor_m = params["proposal_factor_m"];
  const double proposal_factor_C = params["proposal_factor_C"];
  const int locs_dim = locs.n_cols;
  const int p = 2*J*J+1;
  const int T = Y.n_cols;
  const int S = Y.n_rows;

  // Create matrices and cubes for FFBS
  Y.insert_cols(0, 1); // make Y true-indexed; i.e. index 1 is t_1
  arma::cube theta(p, T+1, n_samples), G(p, p, n_samples+1);
  theta.slice(0).zeros();
  arma::mat a(p, T+1), m(p, T+1);
  arma::cube R_inv(p, p, T+1), C(p, p, T+1), C_T(p, p, n_samples+1);
  m.col(0) = m_0;
  C.slice(0) = C_0;

  // Create observation matrix (F) and initial process matrix (G)
  arma::mat w_for_B = makeW(locs, J, L);
  arma::mat F = makeF(locs, w_for_B, J, L);
  const arma::mat FtFiFt = arma::solve(F.t() * F, F.t());
  arma::mat B(S, 2*J*J + 1);
  makeB(B, m_kernel, C_kernel, locs, w_for_B, J, L);
  G.slice(0) = FtFiFt * B;
  
  // Create objects for storing sampled mu_kernel and Sigma_kernel
  arma::mat mu_kernel(locs_dim, n_samples+1);
  arma::cube Sigma_kernel(locs_dim, locs_dim, n_samples+1);
  mu_kernel.col(0) = m_kernel;
  Sigma_kernel.slice(0) = C_kernel;
  arma::colvec mu_kernel_proposal;
  arma::mat Sigma_kernel_proposal, G_proposal;
  arma::mat mu_kernel_proposal_var = proposal_factor_C * C_kernel;
  double mh_ratio;

  // Create variance parameters
  arma::vec sigma2, lambda(n_samples+1);
  lambda.at(0) = rigamma(alpha_lambda, beta_lambda);
  if (sample_sigma2) {
    sigma2.set_size(n_samples+1);
    //sigma2[0] = rigamma(alpha_sigma2, beta_sigma2);
    SampleSigma2(sigma2.at(0), alpha_sigma2, beta_sigma2, Y, F, theta.slice(0));
  }

  // Begin MCMC
  for (int i=0; i<n_samples; ++i) {
    checkUserInterrupt();

    // FFBS
    if (verbose) {
      Rcout << "Filtering sample number " << i+1 << std::endl;
    }
    if (sample_sigma2) sigma2_i = sigma2.at(i);
    KalmanDiscount(m, C, a, R_inv, Y, F, G.slice(i), sigma2_i, lambda.at(i));


    if (verbose) {
      Rcout << "Drawing sample number " << i+1 << std::endl;
    }
    BackwardSample(theta.slice(i), m, a, C, G.slice(i), R_inv);
    C_T.slice(i+1) = C.slice(T); // Save for predictions

    // Sigma2
    if (sample_sigma2) {
      SampleSigma2(sigma2.at(i+1), alpha_sigma2, beta_sigma2, Y, F, theta.slice(i));
    }

    // Lambda (W)
    SampleLambda(lambda.at(i+1), alpha_lambda, beta_lambda,
                 G.slice(i), C, theta.slice(i));

    // MH step for mu
    // Sample proposal value
    mu_kernel_proposal = mvnorm(mu_kernel.col(i), mu_kernel_proposal_var);
    makeB(B, mu_kernel_proposal, Sigma_kernel.slice(i), locs, w_for_B, J, L);
    G_proposal = FtFiFt * B; 
    mh_ratio = 0;
    mh_ratio += ldmvnorm(mu_kernel_proposal, m_kernel, C_kernel);
    mh_ratio -= ldmvnorm(mu_kernel.col(i), m_kernel, C_kernel);
    mh_ratio += kernelLikelihood(G_proposal, theta.slice(i), C);
    mh_ratio -= kernelLikelihood(G.slice(i), theta.slice(i), C);
    mh_ratio -= ldmvnorm(mu_kernel_proposal, mu_kernel.col(i), mu_kernel_proposal_var);
    mh_ratio += ldmvnorm(mu_kernel.col(i), mu_kernel_proposal, mu_kernel_proposal_var);
    if ( std::log(R::runif(0, 1) ) < mh_ratio) {
      mu_kernel.col(i+1) = mu_kernel_proposal;
      G.slice(i+1) = G_proposal;
    } else {
      G.slice(i+1) = G.slice(i);
    }
    
    // Create new B and G matrix from value
    // Likelihood looks something like this
    // (theta_t - G * theta_{t-1})' (lambda * C)^{-1} (theta_t - G * theta_{t-1})
    // Where C comes from the filtering recursions

    // MH step for Sigma
    // Similar to portion for mu

    // Make new process matrix
    //arma::mat B = makeB(m_0, C_0, locs, w_for_B, J, L);
    G.slice(i+1) = FtFiFt * B;
  }

  List results;
  results["m"] = m;
  results["C"] = C;
  results["a"] = a;
  results["theta"]  = theta;
  results["lambda"] = lambda;
  if (sample_sigma2) {
    results["sigma2"] = sigma2;
  } else {
    results["sigma2"] = sigma2_i;
  }
  results["G"] = G;
  results["F"] = F;
  results["C_T"] = C_T;
  return results;
}
