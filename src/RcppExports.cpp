// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// eof
List eof(arma::mat Y, arma::mat F, arma::mat G_0, arma::mat Sigma_G_inv, arma::colvec m_0, arma::mat C_0, arma::mat scale_W, NumericVector params, CharacterVector proc_model, const int n_samples, const bool verbose);
RcppExport SEXP _ideq_eof(SEXP YSEXP, SEXP FSEXP, SEXP G_0SEXP, SEXP Sigma_G_invSEXP, SEXP m_0SEXP, SEXP C_0SEXP, SEXP scale_WSEXP, SEXP paramsSEXP, SEXP proc_modelSEXP, SEXP n_samplesSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G_0(G_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_G_inv(Sigma_G_invSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type m_0(m_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C_0(C_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type scale_W(scale_WSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type proc_model(proc_modelSEXP);
    Rcpp::traits::input_parameter< const int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(eof(Y, F, G_0, Sigma_G_inv, m_0, C_0, scale_W, params, proc_model, n_samples, verbose));
    return rcpp_result_gen;
END_RCPP
}
// ide
List ide(arma::mat Y, arma::mat locs, arma::colvec m_0, arma::mat C_0, arma::mat mean_mu_kernel, arma::mat var_mu_kernel, arma::cube K, arma::cube scale_Sigma_kernel, arma::mat scale_W, NumericVector params, const int n_samples, const bool verbose);
RcppExport SEXP _ideq_ide(SEXP YSEXP, SEXP locsSEXP, SEXP m_0SEXP, SEXP C_0SEXP, SEXP mean_mu_kernelSEXP, SEXP var_mu_kernelSEXP, SEXP KSEXP, SEXP scale_Sigma_kernelSEXP, SEXP scale_WSEXP, SEXP paramsSEXP, SEXP n_samplesSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type m_0(m_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C_0(C_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mean_mu_kernel(mean_mu_kernelSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var_mu_kernel(var_mu_kernelSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type scale_Sigma_kernel(scale_Sigma_kernelSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type scale_W(scale_WSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(ide(Y, locs, m_0, C_0, mean_mu_kernel, var_mu_kernel, K, scale_Sigma_kernel, scale_W, params, n_samples, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ideq_eof", (DL_FUNC) &_ideq_eof, 11},
    {"_ideq_ide", (DL_FUNC) &_ideq_ide, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_ideq(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
