#ifndef KALMAN_H
#define KALMAN_H

#include <RcppArmadillo.h>
using namespace Rcpp;

void CheckDims(arma::mat & Y, arma::mat & F, arma::mat & V,
                 arma::mat & G, arma::mat & W,
                 arma::colvec & m_0, arma::mat & C_0,
                 const int & T, const int & j, const int & p);

arma::colvec mvnorm(arma::colvec M, arma::mat C);

void Kalman(arma::mat & Y, arma::mat & F, arma::mat & V,
                 arma::mat & G, arma::mat & W,
                 arma::mat & m, arma::cube & C,
                 const int & T, const int & j,
                 arma::mat & a, arma::cube & R);



#endif
