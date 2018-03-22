#ifndef KALMAN_H
#define KALMAN_H

#include <RcppArmadillo.h>
using namespace Rcpp;

void CheckDims(arma::mat & Y, arma::mat & F, arma::mat & V,
                 arma::mat & G, arma::mat & W,
                 arma::colvec & m_0, arma::mat & C_0,
                 const int & T, const int & S, const int & p);

void Kalman(arma::mat & Y, arma::mat & F, arma::mat & V,
                 arma::mat & G, arma::mat & W,
                 arma::mat & m, arma::cube & C,
                 const int & T, const int & S,
                 arma::mat & a, arma::cube & R);

void KalmanDiscounted(arma::mat & Y, arma::mat & F, arma::mat & G,
                      arma::mat & m, arma::cube & C,
                      arma::mat & a, arma::cube & R,
                      const int T, const int S, const int p,
                      double sigma2 , double lambda);

#endif
