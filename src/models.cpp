#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <rgen.h>
// [[Rcpp::depends(rgen)]]

#include "kalman.h"
#include "sample.h"
#include "distributions.h"
#include "ide_helpers.h"
#include "misc_helpers.h"

using namespace Rcpp;

//' Fits a DSTM using a wishart prior for W
//'
//' @export
//' @examples
//' # Duhh...nothing yet
//' @importFrom Rcpp sourceCpp evalCpp
//' @useDynLib ideq
// [[Rcpp::export]]
List eof(arma::mat Y, arma::mat F, arma::mat G_0, arma::mat Sigma_G_inv,
         arma::colvec m_0, arma::mat C_0, arma::mat C_W,
         NumericVector params, CharacterVector proc_model,
         const int n_samples, const bool verbose) {
  
  // Extract scalar parameters
  bool AR = proc_model(0) == "AR";
  bool DENSE = proc_model(0) == "Dense";
  const int P = G_0.n_rows;
  const int T = Y.n_cols;
  const int S = Y.n_rows;
  const double alpha_sigma2 = params["alpha_sigma2"];
  const double beta_sigma2  = params["beta_sigma2"];
  double sigma2_i = params["sigma2"];
  const bool sample_sigma2  = sigma2_i < 0;
  const double alpha_lambda = params["alpha_lambda"];
  const double beta_lambda  = params["beta_lambda"];
  const double df_W = params["df_W"];
  const bool discount = C_W.at(0, 0) == NA;
  
  // Create matrices and cubes for FFBS
  Y.insert_cols(0, 1); // make Y true-indexed; i.e. index 1 is t_1
  arma::cube theta(P, T+1, n_samples), G;
  theta.slice(0).zeros();
  arma::mat a(P, T+1), m(P, T+1), tmp;
  arma::cube R_inv(P, P, T+1), C(P, P, T+1), C_T; // C_T for discount models only
  m.col(0) = m_0;
  C.slice(0) = C_0;
  
  // Process model
  if (AR) {
    G.set_size(P, P, n_samples+1);
    G.zeros();
    tmp = G_0.diag();
    G_0.set_size(P, 1);
    G_0 = tmp;
    G.slice(0).diag() = rmvnorm(G_0, arma::inv_sympd(Sigma_G_inv));
  } else if (DENSE) {
    G.set_size(P, P, n_samples+1);
    G_0.reshape(P*P, 1);
    tmp = rmvnorm(G_0, arma::inv_sympd(Sigma_G_inv));
    tmp.reshape(P, P);
    G.slice(0) = tmp;
  } else {
    G.set_size(P, P, 1);
    G.slice(0) = G_0;
  }
  
  // Observation error
  arma::vec sigma2;
  if (sample_sigma2) {
    sigma2.set_size(n_samples+1);
    sampleSigma2(sigma2.at(0), alpha_sigma2, beta_sigma2, Y, F, theta.slice(0));
  }
  else {
    sigma2_i = params["sigma2"];
  }
  
  // Process error
  arma::vec lambda;
  arma::cube W;
  if (discount) {
    lambda.set_size(n_samples+1);
    lambda.at(0) = rigamma(alpha_lambda, beta_lambda);
    C_T.set_size(P, P, n_samples+1);
  } 
  else {
    W.set_size(P, P, n_samples+1);
    W.slice(0) = df_W * C_W;
  }
  
  // Sampling loop
  int G_idx = 0; // This value is incremented each iteration for AR and Dense models
  for (int i = 0; i < n_samples; ++i) {
    checkUserInterrupt();
    
    // Set sigma2_i for FFBS
    if (sample_sigma2) sigma2_i = sigma2.at(i);
    
    // FFBS
    if (verbose) {
      Rcout << "Filtering sample number " << i+1 << std::endl;
    }
    if (discount) {
      kalmanDiscount(m, C, a, R_inv, Y, F, G.slice(G_idx), sigma2_i, lambda.at(i));
      C_T.slice(i+1) = C.slice(T); // Save for predictions
    } else {
      kalman(m, C, a, R_inv, Y, F, G.slice(G_idx), sigma2_i, W.slice(i));
    }
    
    if (verbose) {
      Rcout << "Drawing sample number " << i+1 << std::endl;
    }
    backwardSample(theta.slice(i), m, a, C, G.slice(G_idx), R_inv);
    
    // Sample Sigma2
    if (sample_sigma2) {
      sampleSigma2(sigma2.at(i+1), alpha_sigma2, beta_sigma2, Y, F, theta.slice(i));
    }
    
    // Sample W
    if (discount) {
      sampleLambda(lambda.at(i+1), alpha_lambda, beta_lambda,
                   G.slice(G_idx), C, theta.slice(i));
    } else {
      sampleW(W.slice(i+1), theta.slice(i), G.slice(G_idx), C_W, df_W);
    }
    
    // Sample G
    if (AR) {
      if (discount) {
        sampleAR(G.slice(G_idx+1), R_inv, theta.slice(i),
                 Sigma_G_inv, G_0, true, lambda.at(i+1));
      } else {
        sampleAR(G.slice(G_idx+1), W.slices(i+1, i+1), theta.slice(i), Sigma_G_inv, G_0);
      }
    } else if (DENSE) {
      if (discount) {
        sampleG(G.slice(G_idx+1), R_inv, theta.slice(i),
                Sigma_G_inv, G_0, true, lambda.at(i+1));
      } else {
        sampleG(G.slice(G_idx+1), W.slices(i+1, i+1), theta.slice(i), Sigma_G_inv, G_0);
      }
    }
    
    // Increment G_idx for AR and Dense models
    if (AR || DENSE) ++G_idx;
  }
  
  List results;
  results["F"] = F;
  results["theta"]  = theta;
  if (sample_sigma2) {    results["sigma2"] = sigma2;
  } else {
    results["sigma2"] = sigma2_i;
  }
  if (AR || DENSE) {
    results["G"] = G;
  }
  if (discount) {
    results["lambda"] = lambda;
  } else {
    results["W"] = W;
  }
  
  return results;
};

//' Fits an integrodifference equation model (IDE)
//'
//' @export
//' @examples @importFrom Rcpp sour
//' # Duhh...nothing yet
//'ceCpp evalCpp
//' @useDynLib ideq
// [[Rcpp::export]]
List ide(arma::mat Y, arma::mat locs, arma::colvec m_0, arma::mat C_0,
         arma::colvec mu_kernel_mean, arma::mat mu_kernel_var, arma::cube K,
         arma::cube Sigma_kernel_scale, arma::mat C_W, NumericVector params, 
         const int n_samples, const bool verbose) {
  // Extract scalar parameters
  const double J  = params["J"];
  const double L  = params["L"];
  const int P = 2*J*J+1;
  const int T = Y.n_cols;
  const int S = Y.n_rows;
  const int locs_dim = locs.n_cols;
  const int n_knots = K.n_cols;
  const double alpha_sigma2  = params["alpha_sigma2"];
  const double beta_sigma2   = params["beta_sigma2"];
  double       sigma2_i      = params["sigma2"];
  const bool   sample_sigma2 = sigma2_i < 0;
  const double alpha_lambda  = params["alpha_lambda"];
  const double beta_lambda   = params["beta_lambda"];
  const double df_W = params["df_W"];
  const bool discount = C_W.at(0, 0) == NA;
  const float proposal_factor_mu = params["proposal_factor_mu"];
  const double proposal_factor_Sigma = params["proposal_factor_Sigma"];
  const double Sigma_kernel_df = params["Sigma_kernel_df"];
  const bool SV = params["SV"] > 0;
  const bool dyanamic_K = K.n_slices > 1;
  int K_idx = 0;
  int mu_acceptances = 0;
  int Sigma_acceptances = 0;
  
  // Create matrices and cubes for FFBS
  Y.insert_cols(0, 1); // make Y true-indexed; i.e. index 1 is t_1
  arma::cube theta(P, T+1, n_samples), G(P, P, n_samples+1);
  theta.slice(0).zeros();
  arma::mat a(P, T+1), m(P, T+1), tmp;
  arma::cube R_inv(P, P, T+1), C(P, P, T+1), C_T(P, P, n_samples+1);
  m.col(0) = m_0;
  C.slice(0) = C_0;
  
  // Create objects for storing sampled mu_kernel and Sigma_kernel
  float u;
  arma::cube mu_kernel, mu_kernel_knots;
  arma::cube Sigma_kernel_proposal, Sigma_kernel_knots_proposal;
  arma::field<arma::cube> Sigma_kernel(n_samples+1);
  arma::field<arma::cube> Sigma_kernel_knots(n_samples+1);
  
  if (SV) {
    // Set size of kernel parameter objects
    mu_kernel.set_size(S, locs_dim, n_samples+1);
    mu_kernel_knots.set_size(n_knots, locs_dim, n_samples+1);
    Sigma_kernel_proposal.set_size(locs_dim, locs_dim, S);
    Sigma_kernel_knots_proposal.set_size(locs_dim, locs_dim, n_knots);
    for (int i=0; i<=n_samples; ++i) {
      Sigma_kernel.at(i).set_size(locs_dim, locs_dim, S);
      Sigma_kernel_knots.at(i).set_size(locs_dim, locs_dim, n_knots);
    }
    
    // Set initial values
    mu_kernel_knots.slice(0) = arma::reshape(mu_kernel_mean, n_knots, locs_dim);
    mu_kernel.slice(0) = K.slice(0) * mu_kernel_knots.slice(0);
    Sigma_kernel_knots.at(0) = Sigma_kernel_scale / (Sigma_kernel_df-locs_dim-1);
    mapSigma(Sigma_kernel.at(0), Sigma_kernel_knots.at(0), K);
  }
  else {
    mu_kernel.set_size(locs_dim, 1, n_samples+1);
    mu_kernel.slice(0) = mu_kernel_mean;
    Sigma_kernel_proposal.set_size(locs_dim, locs_dim, 1);
    for (int i=0; i<=n_samples; ++i) {
      Sigma_kernel.at(i).set_size(locs_dim, locs_dim, 1);
    }
    Sigma_kernel.at(0) = Sigma_kernel_scale / (Sigma_kernel_df-locs_dim-1);
  }
  
  arma::mat mu_kernel_proposal, mu_kernel_knots_proposal, G_proposal;
  arma::mat mu_kernel_proposal_var = proposal_factor_mu * proposal_factor_mu
                                                        * mu_kernel_var;
  double Sigma_kernel_proposal_df = locs_dim + Sigma_kernel_df/proposal_factor_Sigma;
  const double Sigma_kernel_adjustment = Sigma_kernel_proposal_df - locs_dim - 1;
  double mh_ratio;
  
  // Create observation matrix (F) and initial process matrix (G)
  arma::mat w_for_B = makeW(J, L);
  arma::mat F = makeF(locs, w_for_B, J, L);
  const arma::mat FtFiFt = arma::solve(F.t() * F, F.t());
  arma::mat B(S, 2*J*J + 1);
  makeB(B, mu_kernel.slice(0), Sigma_kernel.at(0), locs, w_for_B, J, L);
  G.slice(0) = FtFiFt * B;
  
  // Observation error
  arma::vec sigma2;
  if (sample_sigma2) {
    sigma2.set_size(n_samples+1);
    sampleSigma2(sigma2.at(0), alpha_sigma2, beta_sigma2, Y, F, theta.slice(0));
  }
  
  // Process error
  arma::vec lambda;
  arma::cube W;
  if (discount) {
    lambda.set_size(n_samples+1);
    lambda.at(0) = rigamma(alpha_lambda, beta_lambda);
  } 
  else {
    W.set_size(P, P, n_samples+1);
    W.slice(0) = df_W * C_W;
  }
  
  // Begin MCMC
  for (int i=0; i<n_samples; ++i) {
    checkUserInterrupt();
    
    // FFBS
    if (verbose) {
      Rcout << "Filtering sample number " << i+1 << std::endl;
    }
    if (sample_sigma2) sigma2_i = sigma2.at(i);
    
    if (discount) {
      kalmanDiscount(m, C, a, R_inv, Y, F, G.slice(i), sigma2_i, lambda.at(i));
    } else {
      kalman(m, C, a, R_inv, Y, F, G.slice(i), sigma2_i, W.slice(i));
    }
    
    if (verbose) {
      Rcout << "Drawing sample number " << i+1 << std::endl;
    }
    backwardSample(theta.slice(i), m, a, C, G.slice(i), R_inv);
    C_T.slice(i+1) = C.slice(T); // Save for predictions
    
    // Sigma2
    if (sample_sigma2) {
      sampleSigma2(sigma2.at(i+1), alpha_sigma2, beta_sigma2, Y, F, theta.slice(i));
    }
    
    // Process error
    if (discount) {
      sampleLambda(lambda.at(i+1), alpha_lambda, beta_lambda,
                   G.slice(i), C, theta.slice(i));}
    else {
      sampleW(W.slice(i+1), theta.slice(i), G.slice(i), C_W, df_W);
    }
    
    // MH step for mu
    // Propose mu and calculate probabilities under prior
    if (SV) {
      mu_kernel_knots_proposal = proposeMu(mu_kernel_knots.slice(i), mu_kernel_proposal_var);
      mh_ratio  = ldmvnorm(mu_kernel_knots_proposal, mu_kernel_knots.slice(0),
                           mu_kernel_var);
      mh_ratio -= ldmvnorm(mu_kernel_knots.slice(i), mu_kernel_knots.slice(0),
                           mu_kernel_var);
      mu_kernel_proposal = K.slice(0) * mu_kernel_knots_proposal; // map to spatial locations
    } else {
      mu_kernel_proposal = proposeMu(mu_kernel.slice(i), mu_kernel_proposal_var);
      mh_ratio  = ldmvnorm(mu_kernel_proposal, mu_kernel.slice(0), mu_kernel_var);
      mh_ratio -= ldmvnorm(mu_kernel.slice(i), mu_kernel.slice(0), mu_kernel_var);
    }
    
    // Calculate implied proposal for G and likelihood
    makeB(B, mu_kernel_proposal, Sigma_kernel.at(i), locs, w_for_B, J, L);
    G_proposal = FtFiFt * B; 
    
    if (discount) {
      mh_ratio += kernelLikelihoodDiscount(G_proposal, theta.slice(i), C, lambda.at(i+1));
      mh_ratio -= kernelLikelihoodDiscount(G.slice(i), theta.slice(i), C, lambda.at(i+1));
    } else {
      mh_ratio += kernelLikelihood(G_proposal, theta.slice(i), W.slice(i+1));
      mh_ratio -= kernelLikelihood(G.slice(i), theta.slice(i), W.slice(i+1));
    }
    
    // Accept according to mh-ratio
    u = R::runif(0, 1);
    if (std::log(u) < mh_ratio) {
      ++mu_acceptances;
      mu_kernel.slice(i+1) = mu_kernel_proposal;
      G.slice(i+1) = G_proposal;
      if (SV) mu_kernel_knots.slice(i+1) = mu_kernel_knots_proposal;
    } else {
      mu_kernel.slice(i+1) = mu_kernel.slice(i);
      G.slice(i+1) = G.slice(i);
      if (SV) mu_kernel_knots.slice(i+1) = mu_kernel_knots.slice(i);
    }
    
    // MH step for Sigma
    // Propose Sigma and calculate probabilities under prior
    if (SV) {
      for (int k=0; k<n_knots; ++k) {
        Sigma_kernel_knots_proposal.slice(k) = rgen::riwishart(
                                                 Sigma_kernel_proposal_df,
                                                 Sigma_kernel_knots.at(i).slice(k) *
                                                 Sigma_kernel_adjustment);
      }
      mh_ratio  = ldiwishart(Sigma_kernel_knots_proposal, Sigma_kernel_df, 
                             Sigma_kernel_scale);
      mh_ratio -= ldiwishart(Sigma_kernel_knots.at(i), Sigma_kernel_df,
                             Sigma_kernel_scale);
      mapSigma(Sigma_kernel_proposal, Sigma_kernel_knots_proposal, K.slice(0));
    } else {
      Sigma_kernel_proposal.slice(0) = rgen::riwishart(Sigma_kernel_proposal_df,
                                                       Sigma_kernel.at(i).slice(0) * 
                                                       Sigma_kernel_adjustment);
      mh_ratio  = ldiwishart(Sigma_kernel_proposal, Sigma_kernel_df, 
                             Sigma_kernel_scale);
      mh_ratio -= ldiwishart(Sigma_kernel.at(i), Sigma_kernel_df,
                             Sigma_kernel_scale);
    }
    
    // Calculate likelihood
    makeB(B, mu_kernel.slice(i+1), Sigma_kernel_proposal, locs, w_for_B, J, L);
    G_proposal = FtFiFt * B; 
    
    if (discount) {
      mh_ratio += kernelLikelihoodDiscount(G_proposal, theta.slice(i), C, lambda.at(i+1));
      mh_ratio -= kernelLikelihoodDiscount(G.slice(i), theta.slice(i), C, lambda.at(i+1));
    } else {
      mh_ratio += kernelLikelihood(G_proposal, theta.slice(i), W.slice(i+1));
      mh_ratio -= kernelLikelihood(G.slice(i), theta.slice(i), W.slice(i+1));
    }
    
    // Calculate transition probabilities
    mh_ratio -= ldiwishart(Sigma_kernel_proposal, Sigma_kernel_proposal_df,
                           Sigma_kernel.at(i) * Sigma_kernel_adjustment);
    mh_ratio += ldiwishart(Sigma_kernel.at(i), Sigma_kernel_proposal_df,
                           Sigma_kernel_proposal * Sigma_kernel_adjustment);
    
    // Accept according to mh-ratio
    u = R::runif(0, 1);
    if (std::log(u) < mh_ratio) {
      ++Sigma_acceptances;
      Sigma_kernel.at(i+1) = Sigma_kernel_proposal;
      G.slice(i+1) = G_proposal;
      Sigma_kernel_knots.at(i+1) = Sigma_kernel_knots_proposal; 
    } else {
      Sigma_kernel.at(i+1) = Sigma_kernel.at(i);
      Sigma_kernel_knots.at(i+1) = Sigma_kernel_knots.at(i);
    }
  }
  
  // Drop starting values
  G.shed_slice(0);
  C_T.shed_slice(0);
  mu_kernel.shed_slice(0);
  // Last element of Sigma_kernel is removed in dstm_ide in dstm.R
  
  // Save results to list
  List results;
  results["theta"]  = theta;
  results["G"] = G;
  results["F"] = F;
  results["C_T"] = C_T;
  results["mu_kernel"] = mu_kernel;
  results["Sigma_kernel"] = Sigma_kernel;
  if (discount) {
    results["lambda"] = lambda;
  }
  else {
    W.shed_slice(0);
    results["W"] = W;
  }
  if (sample_sigma2) {
    results["sigma2"] = sigma2;
  } 
  else {
    results["sigma2"] = sigma2_i;
  }
  results["mu_acceptances"] = mu_acceptances;
  results["Sigma_acceptances"] = Sigma_acceptances;
  return results;
};