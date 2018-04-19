#ifndef MODELS_H
#define MODELS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

List dstm_discount(arma::mat Y, arma::mat F, arma::mat G_0, arma::mat Sigma_G_inv,
                   arma::colvec m_0, arma::mat C_0, NumericVector params,
                   CharacterVector proc_model, const int n_samples,
                   const bool verbose);

List dstm_IW(void);

List dstm_IDE(void);

#endif
