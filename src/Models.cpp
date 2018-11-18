#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <rgen.h>
// [[Rcpp::depends(rgen)]]

#include "Kalman.h"
#include "Sample.h"
#include "Distributions.h"
#include "IDE_helpers.h"

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
  bool FULL = proc_model(0) == "Full";
  const int p = G_0.n_rows;
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
  arma::cube theta(p, T+1, n_samples), G;
  theta.slice(0).zeros();
  arma::mat a(p, T+1), m(p, T+1), tmp;
  arma::cube R_inv(p, p, T+1), C(p, p, T+1), C_T; // C_T for discount models only
  m.col(0) = m_0;
  C.slice(0) = C_0;

  // Process model
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

  // Observation error
  arma::vec sigma2;
  if (sample_sigma2) {
    sigma2.set_size(n_samples+1);
    SampleSigma2(sigma2.at(0), alpha_sigma2, beta_sigma2, Y, F, theta.slice(0));
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
    C_T.set_size(p, p, n_samples+1);
  } 
  else {
    W.set_size(p, p, n_samples+1);
    W.slice(0) = df_W * C_W;
  }
  
  // Sampling loop
  int G_idx = 0; // This value is incremented each iteration for AR and Full models
  for (int i = 0; i < n_samples; ++i) {
    checkUserInterrupt();

    // Set sigma2_i for FFBS
   if (sample_sigma2) sigma2_i = sigma2.at(i);
    
    // FFBS
    if (verbose) {
      Rcout << "Filtering sample number " << i+1 << std::endl;
    }
    if (discount) {
      KalmanDiscount(m, C, a, R_inv, Y, F, G.slice(G_idx), sigma2_i, lambda.at(i));
      C_T.slice(i+1) = C.slice(T); // Save for predictions
    } else {
      Kalman(m, C, a, R_inv, Y, F, G.slice(G_idx), sigma2_i, W.slice(i));
    }

    if (verbose) {
      Rcout << "Drawing sample number " << i+1 << std::endl;
    }
    BackwardSample(theta.slice(i), m, a, C, G.slice(G_idx), R_inv);

    // Sample Sigma2
    if (sample_sigma2) {
      SampleSigma2(sigma2.at(i+1), alpha_sigma2, beta_sigma2, Y, F, theta.slice(i));
    }

    // Sample W
    if (discount) {
      SampleLambda(lambda.at(i+1), alpha_lambda, beta_lambda,
                   G.slice(G_idx), C, theta.slice(i));
    } else {
      SampleW(W.slice(i+1), theta.slice(i), G.slice(G_idx), C_W, df_W);
    }

    // Sample G
    if (AR) {
      if (discount) {
        SampleAR(G.slice(G_idx+1), R_inv, theta.slice(i),
                 Sigma_G_inv, G_0, true, lambda.at(i+1));
      } else {
        SampleAR(G.slice(G_idx+1), W.slices(i+1, i+1), theta.slice(i), Sigma_G_inv, G_0);
      }
    } else if (FULL) {
      if (discount) {
        SampleG(G.slice(G_idx+1), R_inv, theta.slice(i),
                Sigma_G_inv, G_0, true, lambda.at(i+1));
      } else {
        SampleG(G.slice(G_idx+1), W.slices(i+1, i+1), theta.slice(i), Sigma_G_inv, G_0);
      }
    }
    ++G_idx;
  }
  
  List results;
  results["F"] = F;
  results["theta"]  = theta;
  if (sample_sigma2) {    results["sigma2"] = sigma2;
  } else {
    results["sigma2"] = sigma2_i;
  }
  if (AR || FULL) {
    results["G"] = G;
  }
  if (discount) {
    results["lambda"] = lambda;
  } else {
    results["W"] = W;
  }
  
  return results;
}

//' Fits an integrodifference equation model (IDE)
//'
//' @export
//' @examples @importFrom Rcpp sour
//' # Duhh...nothing yet
//'ceCpp evalCpp
//' @useDynLib ideq
// [[Rcpp::export]]
List ide_sc(arma::mat Y, arma::mat locs, arma::colvec m_0, arma::mat C_0,
            arma::colvec mu_kernel_mean, arma::mat mu_kernel_var,
            arma::mat Sigma_kernel_scale, arma::mat C_W, NumericVector params, 
            const int n_samples, const bool verbose) {
  // Extract scalar parameters
  const double J  = params["J"];
  const double L  = params["L"];
  const int p = 2*J*J+1;
  const int T = Y.n_cols;
  const int S = Y.n_rows;
  const int locs_dim = locs.n_cols;
  const double alpha_sigma2  = params["alpha_sigma2"];
  const double beta_sigma2   = params["beta_sigma2"];
  double       sigma2_i      = params["sigma2"];
  const bool   sample_sigma2 = sigma2_i < 0;
  const double alpha_lambda  = params["alpha_lambda"];
  const double beta_lambda   = params["beta_lambda"];
  const double df_W = params["df_W"];
  const bool discount = C_W.at(0, 0) == NA;
  const double proposal_factor_mu = params["proposal_factor_mu"];
  const double proposal_factor_Sigma = params["proposal_factor_Sigma"];
  const double Sigma_kernel_df = params["Sigma_kernel_df"];

  // Create matrices and cubes for FFBS
  Y.insert_cols(0, 1); // make Y true-indexed; i.e. index 1 is t_1
  arma::cube theta(p, T+1, n_samples), G(p, p, n_samples+1);
  theta.slice(0).zeros();
  arma::mat a(p, T+1), m(p, T+1);
  arma::cube R_inv(p, p, T+1), C(p, p, T+1), C_T(p, p, n_samples+1);
  m.col(0) = m_0;
  C.slice(0) = C_0;
  
  // Create objects for storing sampled mu_kernel and Sigma_kernel
  arma::mat mu_kernel(locs_dim, n_samples+1);
  arma::cube Sigma_kernel(locs_dim, locs_dim, n_samples+1);
  mu_kernel.col(0) = mu_kernel_mean;
  Sigma_kernel.slice(0) = Sigma_kernel_scale / (Sigma_kernel_df-locs_dim-1);
  arma::colvec mu_kernel_proposal;
  arma::mat Sigma_kernel_proposal, G_proposal;
  arma::mat mu_kernel_proposal_var = sqrt(proposal_factor_mu) * mu_kernel_var;
  double Sigma_kernel_proposal_df = locs_dim + Sigma_kernel_df/proposal_factor_Sigma;
  const double Sigma_kernel_adjustment = Sigma_kernel_proposal_df - locs_dim - 1;
  double mh_ratio;

  // Create observation matrix (F) and initial process matrix (G)
  arma::mat w_for_B = makeW(locs, J, L);
  arma::mat F = makeF(locs, w_for_B, J, L);
  const arma::mat FtFiFt = arma::solve(F.t() * F, F.t());
  arma::mat B(S, 2*J*J + 1);
  makeB(B, mu_kernel_mean, Sigma_kernel.slice(0), locs, w_for_B, J, L);
  G.slice(0) = FtFiFt * B;
  
  // Observation error
  arma::vec sigma2;
  if (sample_sigma2) {
    sigma2.set_size(n_samples+1);
    SampleSigma2(sigma2.at(0), alpha_sigma2, beta_sigma2, Y, F, theta.slice(0));
  }

  // Process error
  arma::vec lambda;
  arma::cube W;
  if (discount) {
    lambda.set_size(n_samples+1);
    lambda.at(0) = rigamma(alpha_lambda, beta_lambda);
  } 
  else {
    W.set_size(p, p, n_samples+1);
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
      KalmanDiscount(m, C, a, R_inv, Y, F, G.slice(i), sigma2_i, lambda.at(i));
    } else {
      Kalman(m, C, a, R_inv, Y, F, G.slice(i), sigma2_i, W.slice(i));
    }
    
    if (verbose) {
      Rcout << "Drawing sample number " << i+1 << std::endl;
    }
    BackwardSample(theta.slice(i), m, a, C, G.slice(i), R_inv);
    C_T.slice(i+1) = C.slice(T); // Save for predictions

    // Sigma2
    if (sample_sigma2) {
      SampleSigma2(sigma2.at(i+1), alpha_sigma2, beta_sigma2, Y, F, theta.slice(i));
    }
    
    // Process error
    if (discount) {
      SampleLambda(lambda.at(i+1), alpha_lambda, beta_lambda,
                 G.slice(i), C, theta.slice(i));}
    else {
      SampleW(W.slice(i+1), theta.slice(i), G.slice(i), C_W, df_W);
    }

    // MH step for mu
    mu_kernel_proposal = mvnorm(mu_kernel.col(i), mu_kernel_proposal_var);
    makeB(B, mu_kernel_proposal, Sigma_kernel.slice(i), locs, w_for_B, J, L);
    G_proposal = FtFiFt * B; 
    mh_ratio  = ldmvnorm(mu_kernel_proposal, mu_kernel_mean, mu_kernel_var);
    mh_ratio -= ldmvnorm(mu_kernel.col(i), mu_kernel_mean, mu_kernel_var);
    
    if (discount) {
      mh_ratio += kernelLikelihoodDiscount(G_proposal, theta.slice(i), C, lambda.at(i+1));
      mh_ratio -= kernelLikelihoodDiscount(G.slice(i), theta.slice(i), C, lambda.at(i+1));
    } else {
      mh_ratio += kernelLikelihood(G_proposal, theta.slice(i), W.slice(i+1));
      mh_ratio -= kernelLikelihood(G.slice(i), theta.slice(i), W.slice(i+1));
    }
    
    if ( log(R::runif(0, 1)) < mh_ratio ) {
      mu_kernel.col(i+1) = mu_kernel_proposal;
      G.slice(i+1) = G_proposal;
    } else {
      mu_kernel.col(i+1) = mu_kernel.col(i);
      G.slice(i+1) = G.slice(i);
    }
    
    // MH step for Sigma
    Sigma_kernel.slice(i+1) = Sigma_kernel.slice(i);
    Sigma_kernel_proposal = rgen::riwishart(Sigma_kernel_proposal_df,
                            Sigma_kernel.slice(i) * Sigma_kernel_adjustment);
    makeB(B, mu_kernel.col(i+1), Sigma_kernel_proposal, locs, w_for_B, J, L);
    G_proposal = FtFiFt * B; 
    mh_ratio  = ldiwishart(Sigma_kernel_proposal, Sigma_kernel_df, 
                           Sigma_kernel_scale);
    mh_ratio -= ldiwishart(Sigma_kernel.slice(i), Sigma_kernel_df,
                           Sigma_kernel_scale);
    
    if (discount) {
      mh_ratio += kernelLikelihoodDiscount(G_proposal, theta.slice(i), C, lambda.at(i+1));
      mh_ratio -= kernelLikelihoodDiscount(G.slice(i), theta.slice(i), C, lambda.at(i+1));
    } else {
      mh_ratio += kernelLikelihood(G_proposal, theta.slice(i), W.slice(i+1));
      mh_ratio -= kernelLikelihood(G.slice(i), theta.slice(i), W.slice(i+1));
    }
    
    mh_ratio -= ldiwishart(Sigma_kernel_proposal, Sigma_kernel_proposal_df,
                           Sigma_kernel.slice(i) * Sigma_kernel_adjustment);
    mh_ratio += ldiwishart(Sigma_kernel.slice(i), Sigma_kernel_proposal_df,
                           Sigma_kernel_proposal * Sigma_kernel_adjustment);
    
    if ( log(R::runif(0, 1)) < mh_ratio ) {
      Sigma_kernel.slice(i+1) = Sigma_kernel_proposal;
      G.slice(i+1) = G_proposal;}
    else {
      Sigma_kernel.slice(i+1) = Sigma_kernel.slice(i);
    }
  }

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
    results["W"] = W;
  }
  if (sample_sigma2) {
    results["sigma2"] = sigma2;
  } 
  else {
    results["sigma2"] = sigma2_i;
  }
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
List ide_sv(arma::mat Y, arma::mat locs, arma::colvec m_0, arma::mat C_0,
            arma::colvec mu_kernel_mean, arma::mat mu_kernel_var,
            arma::mat Sigma_kernel_scale, arma::mat C_W, NumericVector params, 
            const int n_samples, const bool verbose) {
  Rcout << "I don't know what to do yet" << std::endl;
  return 42;
}