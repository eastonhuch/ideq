#ifndef IDE_HELPERS_H
#define IDE_HELPERS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat make_w_for_B(const arma::mat & locs, const int J, const int L);

arma::mat makeF(const arma::mat & locs, const arma::mat & w,
                  const int J, const int L);

arma::mat makeB(arma::colvec mu, arma::mat Sigma, const arma::mat & locs,
                const arma::mat & w, const int J, const int L);

#endif