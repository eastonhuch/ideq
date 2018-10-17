#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat makePhi(const arma::mat & locs, const int J, const int L) {
  arma::colvec freqs = 2*PI/L * arma::regspace(1, J);
  arma::mat w(J*J, 2);
  w.col(0) = arma::repmat(freqs, J, 1);
  w.col(1) = arma::repelem(freqs, J, 1);
  arma::mat Jmat = locs.col(0) * w.col(0).t() +
                    locs.col(1) * w.col(1).t();
  arma::mat Phi(Jmat.n_rows, 2*J*J + 1);
  Phi.col(0).fill(0.5);
  Phi.cols(1, L*L) = arma::cos(Jmat);
  Phi.cols(L*L + 1, 2*L*L) = arma::sin(Jmat);
  Phi /= std::sqrt(L);
  return Phi;
}

arma::mat makeB(arma::colvec mu, arma::mat Sigma,
                const arma::mat & locs, const int J, const int L) {
  arma::colvec freqs = 2*PI/L * arma::regspace(1, J);
  arma::mat w(J*J, 2);
  w.col(0) = arma::repmat(freqs, J, 1);
  w.col(1) = arma::repelem(freqs, J, 1);
  arma::mat Jmat1 = (locs.col(0) + mu(0)) * w.col(0).t() +
                    (locs.col(1) + mu(1)) * w.col(1).t();
  arma::colvec Jvec = Sigma.at(0, 0) * arma::square(w.col(0)) +
                      Sigma.at(1, 1) * arma::square(w.col(1)) +
                      Sigma.at(0, 1) * arma::prod(w, 1);
  arma::mat Jmat2 = arma::kron( arma::ones(locs.n_rows, 1), Jvec.t() );
  Jmat2 = arma::exp(-0.5 * Jmat2);
  arma::mat B(Jmat2.n_rows, 2*L*L + 1);
  B.col(0) = Jmat2.col(0);
  B.cols(1, L*L) = Jmat2 % arma::cos(Jmat1);
  B.cols(L*L + 1, 2*L*L) = Jmat2 % arma::sin(Jmat1);
  B /= std::sqrt(L);
  return B;
};
