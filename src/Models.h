#ifndef MODELS_H
#define MODELS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

List dstm_discount(arma::mat Y, arma::mat F, arma::mat G_0, arma::mat Sigma_G_inv,
                   arma::colvec m_0, arma::mat C_0, NumericVector params,
                   CharacterVector proc_model, const int n_samples,
                   const bool verbose);

List dstm_sample_G(arma::mat & Y, const int n_samples, const int p,
                   const bool verbose, const bool sample_sigma2);

List dstm_AR(arma::mat & Y, const int n_samples, const int p,
             const bool verbose, const bool sample_sigma2);

List dstm_IDE(void);

#endif
