// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

void CheckDims(arma::mat & Y, arma::mat & F, arma::mat & V,
                 arma::mat & G, arma::mat & W,
                 arma::colvec & m_0, arma::mat & C_0,
                 const int & T, const int & j, const int & p) {
  // check size of matrices
  if (F.n_rows != j) {
    Rcerr << "F must be j by p" << std::endl;
  }
  if (V.n_rows !=  j || V.n_cols != j) {
    Rcerr << "V must be q by p" << std::endl;
  }
  if (G.n_rows !=  p || G.n_cols != p) {
    Rcerr << "G must be p by p" << std::endl;
    }
  if (W.n_rows !=  p || W.n_cols != p) {
    Rcerr << "W must be p by p" << std::endl;
  }
  if (m_0.n_elem != p) {
    Rcerr << "F must be j by p" << std::endl;
  }
  if (C_0.n_rows !=  p || C_0.n_cols != p) {
    Rcerr << "C_0 must be p by p" << std::endl;
  }

  return;
};

// [[Rcpp::export]]
arma::colvec mvnorm(arma::colvec M, arma::mat C) {
  int n = M.n_elem;
  arma::colvec z(n);
  for (int i = 0; i < n; i++) {
    z.at(i) = R::rnorm(0.0, 1.0);
  }
  arma::colvec x = M + arma::chol(C).t() * z;
  return x;
};


void Kalman(arma::mat & Y, arma::mat & F, arma::mat & V,
                 arma::mat & G, arma::mat & W,
                 arma::mat & m, arma::cube & C,
                 const int & T, const int & j,
                 arma::mat & a, arma::cube & R) {
  // Don't need to keep
  Rcout << "K0" << std::endl;
  arma::mat Q(j, j);
  arma::colvec f(j);
  Rcout << "K1" << std::endl;

  for (int t = 1; t <= T; t++) {
    // One step ahead predictive distribution of theta
    Rcout << "K2a" << std::endl;
    a.col(t - 1) = G * m.col(t - 1);
    Rcout << "K2b" << std::endl;
    R.slice(t - 1) = G * C.slice(t - 1) * G.t() + W;
    Rcout << "K3" << std::endl;

    // One step ahead predictive distribution of Y_t
    f = F * a.col(t - 1);
    Q = F * R.slice(t - 1) * F.t() + V;
    Rcout << "K4" << std::endl;

    // Filtering distribution of theta
    // NOTE: Y.col(t -1) corresponds to Y_t
    m.col(t) = a.col(t - 1) + R.slice(t - 1) * F.t() *
                          solve(Q, (Y.col(t - 1) - f));
    C.slice(t) = R.slice(t - 1) - R.slice(t - 1) * F.t() *
                                  solve(Q, F * R.slice(t - 1));
    Rcout << "K5" << std::endl;
  }
  return;
};

