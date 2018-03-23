  // [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

void CheckDims(arma::mat & Y, arma::mat & F, arma::mat & G,
               arma::colvec & m_0, arma::mat & C_0,
               const int & T, const int & S, const int & p) {
  // check size of matrices
  if (F.n_rows != S) {
    Rcerr << "F must be S by p" << std::endl;
  }
  if (G.n_rows !=  p || G.n_cols != p) {
    Rcerr << "G must be p by p" << std::endl;
    }
  if (m_0.n_elem != p) {
    Rcerr << "m_0 must have p elements" << std::endl;
  }
  if (C_0.n_rows !=  p || C_0.n_cols != p) {
    Rcerr << "C_0 must be p by p" << std::endl;
  }
  return;
};

void Kalman(arma::mat & Y, arma::mat & F, arma::mat & V,
                 arma::mat & G, arma::mat & W,
                 arma::mat & m, arma::cube & C,
                 const int & T, const int & S,
                 arma::mat & a, arma::cube & R) {
  // Don't need to keep these quantities
  arma::mat Q(S, S);
  arma::colvec f(S);

  for (int t = 1; t <= T; ++t) {
    checkUserInterrupt();
    // One step ahead predictive distribution of theta
    a.col(t) = G * m.col(t - 1);
    R.slice(t) = G * C.slice(t - 1) * G.t() + W;

    // One step ahead predictive distribution of Y_t
    f = F * a.col(t);
    Q = F * R.slice(t) * F.t() + V;

    // Filtering distribution of theta
    m.col(t) = a.col(t) + R.slice(t) * F.t() *
                          solve(Q, (Y.col(t) - f));
    C.slice(t) = R.slice(t) - R.slice(t) * F.t() *
                                  solve(Q, F * R.slice(t));
  }
  return;
};

void KalmanDiscounted(arma::mat & Y, arma::mat & F, arma::mat & G,
                      arma::mat & m, arma::cube & C,
                      arma::mat & a, arma::cube & R,
                      const int T, const int S, const int p,
                      double sigma2 , double lambda) {
  // Don't need to keep these quantities
  arma::mat Q(S, S);
  arma::colvec f(S);

  for (int t = 1; t <= T; ++t) {
    checkUserInterrupt();

    // One step ahead predictive distribution of theta
    a.col(t) = G * m.col(t - 1);
    R.slice(t) = (1 + lambda) * G * C.slice(t - 1) * G.t();

    // One step ahead predictive distribution of Y_t
    f = F * a.col(t);
    Q = F * R.slice(t) * F.t();
    Q.diag() += sigma2;

    // Filtering distribution of theta
    m.col(t) = a.col(t) + R.slice(t) * F.t() *
                          solve(Q, (Y.col(t) - f));
    C.slice(t) = R.slice(t) - R.slice(t) * F.t() *
                              solve(Q, F * R.slice(t));
  }
  return;
};
