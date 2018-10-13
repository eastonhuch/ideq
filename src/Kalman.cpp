  // [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <exception>

using namespace Rcpp;

void Kalman(arma::mat & m, arma::cube & C, arma::mat & a, arma::cube & R_inv,
            const arma::mat & Y, const arma::mat & F, const arma::mat & G,
            const arma::mat & W, const double sigma2) {
  // This assumes that V is (sigma2 * I)
  const int T = Y.n_cols - 1;
  const int S = Y.n_rows;

  // Don't need to keep these quantities
  arma::mat Q(S, S);
  arma::colvec f(S);

  for (int t = 1; t <= T; ++t) {
    checkUserInterrupt();
    // One step ahead predictive distribution of theta
    a.col(t) = G * m.col(t - 1);
    R_inv.slice(t) = G * C.slice(t - 1) * G.t() + W;

    // One step ahead predictive distribution of Y_t
    f = F * a.col(t);
    Q = F * R_inv.slice(t) * F.t(); // + V
    Q.diag() += sigma2;

    // Filtering distribution of theta
    m.col(t) = a.col(t) + R_inv.slice(t) * F.t() *
                          solve(Q, (Y.col(t) - f));
    C.slice(t) = R_inv.slice(t) - R_inv.slice(t) * F.t() *
                                  solve(Q, F * R_inv.slice(t));

    // Invert R for sampling
    R_inv.slice(t) = arma::inv_sympd(R_inv.slice(t));
  }
  return;
};

void KalmanDiscount(arma::mat & m, arma::cube & C, arma::mat & a, arma::cube & R_inv,
                    const arma::mat & Y, const arma::mat & F, const arma::mat & G,
                    const double sigma2 , const double lambda) {
  const int T = Y.n_cols - 1;
  const int S = Y.n_rows;

  // Don't need to keep these
  arma::mat Q(S, S);
  arma::colvec f(S);

  for (int t = 1; t <= T; ++t) {
    checkUserInterrupt();

    // One step ahead predictive distribution of theta
    a.col(t) = G * m.col(t - 1);
    R_inv.slice(t) = (1 + lambda) * G * C.slice(t - 1) * G.t();
    // NOTE: R is inverted at the end of this loop

    // One step ahead predictive distribution of Y_t
    f = F * a.col(t);
    Q = F * R_inv.slice(t) * F.t();
    Q.diag() += sigma2;

    // Filtering distribution of theta
    m.col(t) = a.col(t) + R_inv.slice(t) * F.t() *
                          solve(Q, (Y.col(t) - f));
    C.slice(t) = R_inv.slice(t) - R_inv.slice(t) * F.t() *
                              solve(Q, F * R_inv.slice(t));

    // Invert R for sampling
    try {
      R_inv.slice(t) = arma::inv_sympd(R_inv.slice(t));
    }
    catch (std::runtime_error e)
    {
      Rcout << "Failed to invert R using inv_sympd" << std::endl;
      Rcout << "Forcing symmetry using (R + R')/2 and trying again" << std::endl;
      arma::mat R_sym = ( R_inv.slice(t) + R_inv.slice(t).t() ) / 2;
      R_inv.slice(t) = arma::inv_sympd(R_sym);
    }
  }
  return;
};
