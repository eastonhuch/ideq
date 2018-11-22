#ifndef MODELS_H
#define MODELS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

List eof(arma::mat Y, arma::mat F, arma::mat G_0, arma::mat Sigma_G_inv,
         arma::colvec m_0, arma::mat C_0, arma::mat C_W,
         NumericVector params, CharacterVector proc_model,
         const int n_samples, const bool verbose);

List ide(arma::mat Y, arma::mat locs, arma::colvec m_0, arma::mat C_0,
         arma::colvec mu_kernel_mean, arma::mat mu_kernel_var,
         arma::mat Sigma_kernel_scale, arma::mat C_W, NumericVector params, 
         const int n_samples, const bool verbose);

#endif
