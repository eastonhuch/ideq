#ifndef BACKWARDSAMPLE_H
#define BACKWARDSAMPLE_H

#include <RcppArmadillo.h>
using namespace Rcpp;

void BackwardSample(arma::cube & theta, arma::mat & m, arma::mat & a,
                    arma::cube & C, arma::mat & G, arma::cube & R_inv,
                    const int & n_samples, int & start_slice,
                    const bool & verbose);

void SampleSigma2(const double & alpha_sigma2, const double & beta_sigma2,
                 const int & S, const int & T, int i,
                 arma::mat & Y, arma::mat & F_,
                 arma::cube & theta, arma::colvec & sigma2);

void SampleLambda(const double & alpha_lambda, const double & beta_lambda,
                 const int & p, const int & T, int i,
                 arma::mat & G, arma::cube & C,
                 arma::cube & theta, arma::colvec & lambda);

void SampleG(arma::mat & G, arma::cube & W_inv, arma::mat & theta,
             arma::mat & Sigma_g_inv, arma::mat & mu_g, const int & p,
             const int & T, const bool Discount = false,
             const double lambda = 0.0);

void SampleAR(arma::mat & G, arma::cube & W_inv, arma::mat & theta,
              arma::mat & Sigma_G_inv, arma::mat & mu_G, const int & T,
              const int Discount = false, const double lambda = 0.0);

void SampleV_inv (arma::mat & Y, arma::mat & F, arma::cube & theta,
                  arma::cube & V, arma::mat & C_V, const int & df_V,
                  int & i, const int & T);

void SampleW (arma::mat & theta, arma::mat & G, arma::mat & W,
              arma::mat & C_W, const int & df_W, const int & T);

#endif
