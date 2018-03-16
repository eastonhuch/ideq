// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "Distributions.h"
using namespace Rcpp;

void BackwardSample(arma::cube & theta, arma::mat & m, arma::mat & a,
                    arma::cube & C, arma::mat & G,
                    arma::cube & R, const int & T,
                    const int & n_samples, int & start_slice,
                    const bool & verbose, const int & p) {
  arma::colvec h_t(p);
  arma::mat H_t(p, p);

  // NOTE: n_samples is how many samples we want to draw on this function call
  for (int s = start_slice; s < (start_slice + n_samples); s++) {
    if (verbose) {
      Rcout << "Drawing sample number " << s << std::endl;
    }
    checkUserInterrupt();
    theta.slice(s).col(T) = mvnorm(m.col(T), C.slice(T)); // draw theta_T
    // Draw values for theta_{T-1} down to theta_0
    for (int t = T-1; t >= 0; t--) {
      // Mean and variance of theta_t
      h_t = m.col(t) + C.slice(t) * G.t() *
        solve(R.slice(t), (theta.slice(s).col(t + 1) - a.col(t)));

      H_t = C.slice(t) - C.slice(t) * G.t() *
        solve(R.slice(t), G) * C.slice(t);
      // Draw value for theta_t
      theta.slice(s).col(t) = mvnorm(h_t, H_t);
    }
  }
  return;
}
