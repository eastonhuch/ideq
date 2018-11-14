#ifndef IDE_HELPERS_H
#define IDE_HELPERS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat makeW(const arma::mat & locs, const int J, const int L);

arma::mat makeF(const arma::mat & locs, const arma::mat & w,
                  const int J, const int L);

void makeB(arma::mat & B, const arma::colvec mu, const arma::mat Sigma, 
           const arma::mat & locs, const arma::mat & w, const int J, const int L);

void makeB_SV(arma::mat & B, const arma::mat mu, const arma::cube Sigma, 
              const arma::mat & locs, const arma::mat & w, const int J, const int L);

double kernelLikelihood(const arma::mat & G, const arma::mat & theta, 
                        const arma::cube & C, const double lambda);

#endif
