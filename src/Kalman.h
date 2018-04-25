#ifndef KALMAN_H
#define KALMAN_H

#include <RcppArmadillo.h>
using namespace Rcpp;

void Kalman(const arma::mat & Y, const arma::mat & F, const arma::mat & G,
            const arma::mat & W, arma::mat & m, arma::cube & C, arma::mat & a,
            arma::cube & R_inv, const double sigma2);

void KalmanDiscount(const arma::mat & Y, const arma::mat & F, const arma::mat & G,
                    arma::mat & m, arma::cube & C,
                    arma::mat & a, arma::cube & R_inv,
                    const double sigma2 , const double lambda);

#endif
