#ifndef IDE_HELPERS_H
#define IDE_HELPERS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat makeW(const int J, const int L);

arma::mat makeF(const arma::mat & locs, const arma::mat & w,
                  const int J, const int L);

void makeB(arma::mat & B, const arma::mat & mu, const arma::cube & Sigma, 
           const arma::mat & locs, const arma::mat & w, const int J, const int L);

double kernelLikelihoodDiscount(const arma::mat & G, const arma::mat & theta, 
                        const arma::cube & C, const double lambda);

double kernelLikelihood(const arma::mat & G, const arma::mat & theta, 
                        const arma::mat W);

arma::mat proposeMu(arma::mat mu, arma::mat Sigma);

void mapSigma(arma::cube & s_many, const arma::cube & s_few,
              const arma::mat K);

#endif