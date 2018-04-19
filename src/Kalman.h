#ifndef KALMAN_H
#define KALMAN_H

#include <RcppArmadillo.h>
using namespace Rcpp;

void Kalman(arma::mat & Y, arma::mat & F, arma::mat & G, arma::mat & W,
                 arma::mat & m, arma::cube & C, arma::mat & a,
                 arma::cube & R_inv, const double sigma2);

void KalmanDiscount(arma::mat & Y, arma::mat & F, arma::mat & G,
                      arma::mat & m, arma::cube & C,
                      arma::mat & a, arma::cube & R,
                      const double sigma2 , const double lambda);

#endif
