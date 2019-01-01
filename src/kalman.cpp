#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
  
#include "Misc_helpers.h"

using namespace Rcpp;

void kalman(arma::mat & m, arma::cube & C, arma::mat & a, arma::cube & R_inv,
            const arma::mat & Y, const arma::mat & F, const arma::mat & G,
            const double sigma2, const arma::mat & W) {
  // This assumes that V is (sigma2 * I)
  const int T = Y.n_cols-1;
  const int S = Y.n_rows;
  const int p = G.n_rows;

  // Don't need to keep these quantities
  arma::mat Q(S, S), Q_inv(S, S), R_t(p, p), FR(S, p), RF_t(p, S);
  arma::colvec f(S);
  // Could also save transposed copies of F, G

  for (int t = 1; t <= T; ++t) {
    checkUserInterrupt();
    // One step ahead predictive distribution of theta
    a.col(t) = G * m.col(t-1);
    R_t = G * C.slice(t-1) * G.t() + W;

    // One step ahead predictive distribution of Y_t
    f = F * a.col(t);
    FR = F * R_t;
    Q = FR * F.t();
    Q.diag() += sigma2;
    Q_inv = arma::inv_sympd(Q);

    // Filtering distribution of theta
    RF_t = FR.t();
    m.col(t) = a.col(t) + RF_t * Q_inv * (Y.col(t) - f);
    C.slice(t) = R_t - RF_t * Q_inv * FR;
    //makeSymmetric(C.slice(t));

    // Invert R for sampling
    R_inv.slice(t) = arma::inv_sympd(R_t);
  }
  return;
};

void kalmanDiscount(arma::mat & m, arma::cube & C, arma::mat & a, arma::cube & R_inv,
                    const arma::mat & Y, const arma::mat & F, const arma::mat & G,
                    const double sigma2 , const double lambda) {
  // This assumes that V is (sigma2 * I)
  const int T = Y.n_cols-1;
  const int S = Y.n_rows;
  const int p = G.n_rows;

  // Don't need to keep these quantities
  arma::mat Q(S, S), Q_inv(S, S), R_t(p, p), FR(S, p), RF_t(p, S);
  arma::colvec f(S);
  // Could also save transposed copies of F, G

  for (int t = 1; t <= T; ++t) {
    checkUserInterrupt();
    // One step ahead predictive distribution of theta
    a.col(t) = G * m.col(t-1);
    R_t = (1 + lambda) * G * C.slice(t-1) * G.t();
    makeSymmetric(R_t);

    // One step ahead predictive distribution of Y_t
    f = F * a.col(t);
    FR = F * R_t;
    Q = FR * F.t();
    makeSymmetric(Q);
    Q.diag() += sigma2;
    Q_inv = arma::inv_sympd(Q);
    makeSymmetric(Q_inv);

    // Filtering distribution of theta
    RF_t = FR.t();
    m.col(t) = a.col(t) + RF_t * Q_inv * (Y.col(t) - f);
    C.slice(t) = R_t - RF_t * Q_inv * FR;
    makeSymmetric(C.slice(t));

    // Invert R for sampling
    R_inv.slice(t) = arma::inv_sympd(R_t);
  }
  return;
};