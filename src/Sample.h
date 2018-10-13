#ifndef BACKWARDSAMPLE_H
#define BACKWARDSAMPLE_H

#include <RcppArmadillo.h>
using namespace Rcpp;

void BackwardSample(arma::cube & theta, arma::mat & m, arma::mat & a,
                    arma::cube & C, arma::mat & G, arma::cube & R_inv,
                    const int & n_samples, int & start_slice,
                    const bool & verbose);

void SampleSigma2(double & sigma2_new, const double & alpha_sigma2, const double & beta_sigma2,
                  const arma::mat & Y, const arma::mat & F, const arma::mat & theta);

void SampleLambda(double & lambda_new, const double & alpha_lambda, const double & beta_lambda,
                  const arma::mat & G, const arma::cube & C, const arma::mat & theta);

void SampleG(arma::mat & G, arma::cube & W_inv, arma::mat & theta,
             arma::mat & Sigma_g_inv, arma::mat & mu_g, const int & p,
             const int & T, const bool Discount = false,
             const double lambda = 0.0);

void SampleAR(arma::mat & G, arma::cube & W_inv, arma::mat & theta,
              arma::mat & Sigma_G_inv, arma::mat & mu_G,
              const bool Discount = false, const double lambda = 0.0);

void SampleV(arma::mat & Y, arma::mat & F, arma::cube & theta,
             arma::mat & V, arma::mat & C_V, const int df_V);

void SampleW(arma::mat & theta, arma::mat & G, arma::mat & W,
             arma::mat & C_W, const int df_W);

#endif
