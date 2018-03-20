#ifndef BACKWARDSAMPLE_H
#define BACKWARDSAMPLE_H

#include <RcppArmadillo.h>
using namespace Rcpp;

void BackwardSample(arma::cube & theta, arma::mat & m, arma::mat & a,
                    arma::cube & C, arma::mat & G,
                    arma::cube & R, const int & T,
                    const int & n_samples, int & start_slice,
                    const bool & verbose, const int & p);

void SampleSigma(const double & alpha_sigma, const double & beta_sigma,
                 const int & S, const int & T, int i,
                 arma::mat & Y, arma::mat & F_, arma::mat & a, arma::colvec & sigma);

void SampleTau(const double & alpha_tau, const double & beta_tau,
                 const int & p, const int & T, int i,
                 arma::mat & G, arma::cube & C,
                 arma::mat & a, arma::colvec & tau);
#endif
