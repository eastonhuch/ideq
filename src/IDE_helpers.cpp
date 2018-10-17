#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat make_w_for_B(const arma::mat & locs, const int J, const int L) {
  arma::colvec freqs = 2*PI/L * arma::regspace(1, J);
  arma::mat w(J*J, 2);
  w.col(0) = arma::repmat(freqs, J, 1);
  w.col(1) = arma::repelem(freqs, J, 1);
  return w;
}

// Observation matrix
// The number of total basis function is J^2+1
// L is the range of the Fourier approximation
// locs are the centered/scaled spatial locations
arma::mat makeF(const arma::mat & locs, const arma::mat & w,
                  const int J, const int L) {
  arma::mat Jmat = locs.col(0) * w.col(0).t() +
                   locs.col(1) * w.col(1).t();
  arma::mat Phi(Jmat.n_rows, 2*J*J + 1);
  Phi.col(0).fill(0.5);
  Phi.cols(1, J*J) = arma::cos(Jmat);
  Phi.cols(J*J + 1, 2*J*J) = arma::sin(Jmat);
  Phi /= std::sqrt(L);
  return Phi;
}

// The function makeB returns the matrix B used as part of the process matrix
// mu and Sigma are the parameters of the IDE kernel
arma::mat makeB(arma::colvec mu, arma::mat Sigma, const arma::mat & locs,
                const arma::mat & w, const int J, const int L) {
  arma::mat Jmat1 = (locs.col(0) + mu(0)) * w.col(0).t() +
                    (locs.col(1) + mu(1)) * w.col(1).t();
  arma::colvec Jvec = Sigma.at(0, 0) * arma::square(w.col(0)) +
                      Sigma.at(1, 1) * arma::square(w.col(1)) +
                      Sigma.at(0, 1) * arma::prod(w, 1);
  arma::mat Jmat2 = arma::kron( arma::ones(locs.n_rows, 1), Jvec.t() );
  Jmat2 = arma::exp(-0.5 * Jmat2);
  arma::mat B(Jmat2.n_rows, 2*J*J + 1);
  B.col(0) = Jmat2.col(0);
  B.cols(1, J*J) = Jmat2 % arma::cos(Jmat1);
  B.cols(J*J + 1, 2*J*J) = Jmat2 % arma::sin(Jmat1);
  B /= std::sqrt(L);
  return B;
};
