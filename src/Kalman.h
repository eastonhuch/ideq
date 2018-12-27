#ifndef KALMAN_H
#define KALMAN_H

#include <RcppArmadillo.h>
using namespace Rcpp;

void kalman(arma::mat & m, arma::cube & C, arma::mat & a, arma::cube & R_inv,
            const arma::mat & Y, const arma::mat & F, const arma::mat & G,
            const double sigma2, const arma::mat & W);

void kalmanDiscount(arma::mat & m, arma::cube & C, arma::mat & a, arma::cube & R_inv,
                    const arma::mat & Y, const arma::mat & F, const arma::mat & G,
                    const double sigma2 , const double lambda);

#endif
