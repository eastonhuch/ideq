#ifndef MISC_HELPERS_H
#define MISC_HELPERS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

void make_symmetric(arma::mat & X);

void mapSigma(arma::cube & s_many, const arma::cube & s_few,
              const arma::mat K);

#endif
