#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::colvec mvnorm(const arma::colvec & mean, const arma::mat & Sigma);

double ldmvnorm(const arma::colvec x, const arma::colvec & mean,
                const arma::mat & Sigma);

double ldiwishart(const arma::mat x, const double df,
                  const arma::mat & scale);

double rigamma(const double a, const double scl);

#endif
