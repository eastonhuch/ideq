#' Fits a dynamic spatio-temporal model (DSTM)
#'
#' @param Y S by T matrix containing response variable at S spatial locations and T time points
#' @param model character string; options include `discount`, `sample_G`, `AR`, and `IDE
#' @param n_samples integer; number of posterior samples to take
#' @param p integer; dimension of G in the state equation \eqn{\theta_{t+1} = G \theta_{t}}
#' @param verbose boolean; controls verbosity
#' @param sample_sigma2 whether boolean; to sample \eqn{\sigma^2}
#'
#' @keyword IDE, Kalman, Filter
#' @export
#' @examples
#' # Duhh...nothing yet
dstm <- function(Y, params = NULL, model = "discount", n_samples = 1L, p = 20L,
                 verbose = FALSE, sample_sigma2 = TRUE) {
  results <- "No output...whoops"
  if (model == "discount") {
    # Set m_0
    if ("m_0" %in% names(params)) {
      if(length(params[["m_0"]]) != p) stop("m_0 must have length p")
      m_0 <- params[["m_0"]]
    } else {
      message("No prior was provided for m_0 so I am using a vector of zeros")
      m_0 <- rep(0, p)
    }

    # Set C_0
    if ("C_0" %in% names(params)) {
      if (!is.matrix(params[["C_0"]]) || any(dim(params[["C_0"]]) != c(p, p))){
        stop("C_0 must be a p by p matrix")
      }
      C_0 <- params[["C_0"]]
    } else {
      message("No prior was provided for C_0 so I am using an identity matrix ")
      C_0 <- diag(p)
    }

    # Set F
    if ("F" %in% names(params)) {
      if (!is.matrix(params[["F"]]) ||
          nrow(params[["F"]]) != nrow(Y) ||
          ncol(params[["F"]]) != p) {
        stop("F must be n by p")
      }
      F_ <- params[["F"]]
    } else {
      message("No F matrix was provided so I am using empirical orthogonal functions (eigenvectors)")
      F_ <-eigen(cov(t(Y)))$vectors[, 1:p]
    }

    # Set G
    if ("G" %in% names(params)) {
      if (!is.matrix(params[["G"]]) ||
          any(dim(params[["G"]]) != c(p, p))) {
        stop("G must be a p by p matrix")
      }
      G <- params[["G"]]
    } else {
      message("No G matrix was provided so I am using an identity matrix")
      G <-diag(p)
    }

    # Set prior for lambda
    if ("alpha_lambda" %in% names(params)) {
      alpha_lambda <- params[["alpha_lambda"]]
    } else {
      message("alpha_lambda was not provided so I am using 2.25")
      alpha_lambda <- 2.25
    }
    if ("beta_lambda" %in% names(params)) {
      beta_lambda <- params[["beta_lambda"]]
    } else {
      message("beta_lambda was not provided so I am using 0.0625")
      beta_lambda <- 0.0625
    }

    # Set prior for sigma2 or sigma2 itself if not sampling
    if (sample_sigma2) {
      sigma2 <- NULL
      if ("alpha_sigma2" %in% names(params)) {
        alpha_sigma2 <- params[["alpha_sigma2"]]
      } else {
        v <- mean(apply(Y, 1, var))
        alpha_sigma2 <- 2 + v / 4
        message(paste("alpha_sigma2 was not provided so I am using", alpha_sigma2))
      }
      if ("beta_sigma2" %in% names(params)) {
        beta_sigma2 <- params[["beta_sigma2"]]
      } else {
        if (!exists("v")) v <- mean(apply(Y, 1, var))
        beta_sigma2 <- (alpha_sigma2 - 1) * v
        message(paste("beta_sigma2 was not provided so I am using", beta_sigma2))
      }
      if (alpha_sigma2 <= 0 || beta_sigma2 <= 0) {
        stop("alpha_sigma2 or beta_sigma2 is not positive; specify them manually in params")
      }
    } else {
      alpha_sigma2 <- beta_sigma2 <- NULL
      if ("sigma2" %in% names(params)) {
        sigma2 <- params[["sigma2"]]
      } else {
        if (!exists("v")) v <- mean(apply(Y, 1, var))
        sigma2 <- v
        message(paste("sigma2 was not provided so I am using", sigma2))
      }
    }

    scalar_params <- c(alpha_lambda, beta_lambda, alpha_sigma2, beta_sigma2, sigma2)
    results <- dstm_discount(Y, F_, G, m_0, C_0, scalar_params, n_samples, p, sample_sigma2, verbose)

  } else if (model == "sample_G") {
    results <- dstm_sample_G(Y, n_samples, p, verbose, sample_sigma2)

  } else if (model == "AR") {
    results <- dstm_AR(Y, n_samples, p, verbose, sample_sigma2)

  } else if (model == "IDE") {
    results <- dstm_IDE()

  } else {
    stop("Model type not recognized")
  }

  return(results)
}
