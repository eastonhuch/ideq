#ifndef KALMAN_H
#define KALMAN_H

#include <RcppArmadillo.h>
using namespace Rcpp;

void Kalman(arma::mat & m, arma::cube & C, arma::mat & a, arma::cube & R_inv,
            const arma::mat & Y, const arma::mat & F, const arma::mat & G,
            const arma::mat & W, const double sigma2);

void KalmanDiscount(arma::mat & m, arma::cube & C, arma::mat & a, arma::cube & R_inv,
                    const arma::mat & Y, const arma::mat & F, const arma::mat & G,
                    const double sigma2 , const double lambda);

#endif
