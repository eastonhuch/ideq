#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat makeW(const int J, const int L) {
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
void makeB(arma::mat & B, const arma::mat & mu, const arma::cube & Sigma, 
           const arma::mat & locs, const arma::mat & w, const int J, const int L) {
  
  const bool SV = mu.n_cols > 1;
  
  // Jmat1 and Jmat2
  arma::mat Jmat1, Jmat2;
  
  if (SV) {
    Jmat1 = (locs.col(0) + mu.col(0)) * w.col(0).t() +
            (locs.col(1) + mu.col(1)) * w.col(1).t();
    arma::colvec tmp = Sigma.tube(0, 0);
    Jmat2 = tmp * arma::square(w.col(0).t());
    tmp = Sigma.tube(1, 1);
    Jmat2 += tmp * arma::square(w.col(1).t());
    tmp = Sigma.tube(0, 1);
    Jmat2 += tmp * arma::prod(w.t(), 0);
  } else {
    Jmat1 = (locs.col(0) + mu(0)) * w.col(0).t() +
            (locs.col(1) + mu(1)) * w.col(1).t();
    arma::mat Jvec = Sigma.at(0, 0, 0) * arma::square(w.col(0)) +
                     Sigma.at(1, 1, 0) * arma::square(w.col(1)) +
                     Sigma.at(0, 1, 0) * arma::prod(w, 1);
    Jmat2 = arma::kron( arma::ones(locs.n_rows, 1), Jvec.t() );
  }
  
  // Exponentiate Jmat2
  Jmat2 = arma::exp(-0.5 * Jmat2);
  
  // B
  B.col(0) = Jmat2.col(0);
  B.cols(1, J*J) = Jmat2 % arma::cos(Jmat1);
  B.cols(J*J + 1, 2*J*J) = Jmat2 % arma::sin(Jmat1);
  B /= std::sqrt(L);
  return;
};

double kernelLikelihoodDiscount(const arma::mat & G, const arma::mat & theta, 
                        const arma::cube & C, const double lambda) {
  const int T = theta.n_cols-1;
  const int p = theta.n_rows;
  arma::mat tmp = arma::zeros(1, 1);
  arma::colvec d;
  
  for (int t=1; t<=T; ++t) {
    d = theta.col(t) - G * theta.col(t-1);
    tmp += d.t() * arma::solve(lambda * G * C.slice(t) * G.t(), d);
  }
  
  return -tmp(0)/2; 
};

double kernelLikelihood(const arma::mat & G, const arma::mat & theta,
                        const arma::mat W) {
  const int T = theta.n_cols-1;
  const int p = theta.n_rows;
  arma::mat tmp = arma::zeros(1, 1);
  arma::colvec d;
  
  for (int t=1; t<=T; ++t) {
    d = theta.col(t) - G * theta.col(t-1);
    tmp += d.t() * arma::solve(W, d);
  }
  
  return -tmp(0)/2; 
};