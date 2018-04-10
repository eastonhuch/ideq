#ifndef MODELS_H
#define MODELS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

List dstm_discount(arma::mat & Y, const int n_samples, const int p,
                   const bool sample_sigma2, const bool verbose);

List dstm_sample_G(arma::mat & Y, const int n_samples,
                   const int p, const bool verbose);

List dstm_AR(arma::mat & Y, const int n_samples,
                   const int p, const bool verbose);

List dstm_IDE(void);

#endif
