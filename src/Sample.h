#ifndef BACKWARDSAMPLE_H
#define BACKWARDSAMPLE_H

#include <RcppArmadillo.h>
using namespace Rcpp;

void BackwardSample(arma::cube & theta, arma::mat & m, arma::mat & a,
                    arma::cube & C, arma::mat & G,
                    arma::cube & R, const int & T,
                    const int & n_samples, int & start_slice,
                    const bool & verbose, const int & p);

void SampleSigma2(const double & alpha_sigma2, const double & beta_sigma2,
                 const int & S, const int & T, int i,
                 arma::mat & Y, arma::mat & F_,
                 arma::cube & theta, arma::colvec & sigma2);

void SampleLambda(const double & alpha_lambda, const double & beta_lambda,
                 const int & p, const int & T, int i,
                 arma::mat & G, arma::cube & C,
                 arma::cube & theta, arma::colvec & lambda);
#endif
