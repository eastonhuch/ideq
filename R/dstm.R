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
dstm <- function(Y, obs_model = "EOF", proc_model = "RW",
                 proc_error = "discount", p = 10L,
                 n_samples = 1L, sample_sigma2 = TRUE,
                 verbose = FALSE, params = NULL) {
  results <- "No output...whoops"

  # Observation Model; creates F, m_0, C_0
  if (obs_model == "EOF") {
    F_ <-eigen(cov(t(Y)))$vectors[, 1:p]
    # Set m_0
    m_0 <- NULL
    if ("m_0" %in% names(params)) {
      if(length(params[["m_0"]]) != p) stop("m_0 must have length p")
      m_0 <- params[["m_0"]]
    }
    else {
      message("No prior was provided for m_0 so I am using a vector of zeros")
      m_0 <- rep(0, p)
    }

    # Set C_0
    C_0 <- matrix()
    if ("C_0" %in% names(params)) {
      if (!is.matrix(params[["C_0"]]) || any(dim(params[["C_0"]]) != c(p, p))){
        stop("C_0 must be a p by p matrix")
      }
      C_0 <- params[["C_0"]]
    }
    else {
      message("No prior was provided for C_0 so I am using an identity matrix ")
      C_0 <- diag(p)
    }

    }
  else if (obs_model == "IDE") {
    stop("IDE model not implemented")
  }
  else {
    stop("obs_model is invalid")
  }

  # Process Model; creates G_0 (mu_G) and Sigma_G_inv
  Sigma_G_inv <- matrix()
  if (proc_model == "RW") {
    G_0 <- diag(p)

  }
  else if (proc_model == "AR") {
    if ("mu_G" %in% names(params)) {
      if (!is.diagonal.matrix(params$mu_G) ||
          any(dim(params$mu_G) != c(p, p))) {
        stop("mu_G must be a diagonal p by p matrix")
      } else {
        G_0 <- params$mu_G
      }
    } else {
      message("mu_G was not provided, so I am using an identity matrix")
      G_0 <- diag(p)
    }

    if ("Sigma_G_inv" %in% names(params)) {
      if (!is.positive.definite(params$Sigma_G_inv) ||
          !is.symmetric.matrix(params$Sigma_G_inv)) {
        stop("Sigma_G_inv must be symmetric positive definite matrix")
      } else {
        Sigma_G_inv <- params$Sigma_G_inv
      }
    } else {
      message("Sigma_G_inv was not provided, so I am using an identity matrix")
      Sigma_G_inv <- diag(p)
    }

  }
  else if (proc_model == "Full") {
    if ("mu_G" %in% names(params)) {
      if (!is.matrix(params$mu_G) ||
          any(dim(params$mu_G) != c(p, p))) {
        stop("mu_G must be a p by p matrix")
      } else {
        G_0 <- params$mu_G
      }
    } else {
      message("mu_G was not provided, so I am using an identity matrix")
      G_0 <- diag(p)
    }

    if ("Sigma_G_inv" %in% names(params)) {
      if (!is.positive.definite(params$Sigma_G_inv) ||
          !is.symmetric.matrix(params$Sigma_G_inv) ||
          nrow(params$Sigma_G_inv) != p^2) {
        stop("Sigma_G_inv must be a p^2 by p^2 symmetric positive definite matrix")
      } else {
        Sigma_G_inv <- params$Sigma_G_inv
      }
    } else {
      message("Sigma_G_inv was not provided, so I am using an identity matrix")
      Sigma_G_inv <- diag(p^2)
    }

  }
  else {
    stop("proc_model is invalid")
  }

  # Process Error; creates all necessary params (e.g., alpha_lambda, sigma2)
  if (proc_error == "discount") {

    # Set prior for lambda
    alpha_lambda <- beta_lambda <- NULL
    if ("alpha_lambda" %in% names(params)) {
      alpha_lambda <- params[["alpha_lambda"]]
    }
    else {
      message("alpha_lambda was not provided so I am using 2.25")
      alpha_lambda <- 2.25
    }

    if ("beta_lambda" %in% names(params)) {
      beta_lambda <- params[["beta_lambda"]]
    }
    else {
      message("beta_lambda was not provided so I am using 0.0625")
      beta_lambda <- 0.0625
    }

    # Set prior for sigma2 or sigma2 itself if not sampling
    alpha_sigma2 <- beta_sigma2 <- sigma2 <- NULL
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
    }
    else if ("sigma2" %in% names(params)) {
      sigma2 <- params[["sigma2"]]
    }
    else {
      if (!exists("v")) v <- mean(apply(Y, 1, var))
      sigma2 <- v
      message(paste("sigma2 was not provided so I am using", sigma2))
      }

    # Group scalar params into vector
    scalar_params <- c(alpha_lambda, beta_lambda, alpha_sigma2, beta_sigma2, sigma2)

    # Run the model
    results <- dstm_discount(Y, F_, G_0, Sigma_G_inv, m_0, C_0,
                  scalar_params, proc_model, n_samples, verbose)
  }
  else if (proc_error == "IW") {
    stop("IW process error is not yet implemented")
    results <- dstm_IW()
  }
  else {
    stop("I don't know that type of process error")
  }

  return(results)
}
