#' Fits a dynamic spatio-temporal model (DSTM)
#'
#' @param Y S by T matrix containing response variable at S spatial locations and T time points
#' @param model character string; options include `discount`, `sample_G`, `AR`, and `IDE`
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

  # Observation Error; creates alpha_sigma2, beta_sigma2, sigma2
  # NOTE: We could also draw this from a Wishart distribution...nah
  alpha_sigma2 <- beta_sigma2 <- sigma2 <- NULL
  if (sample_sigma2) {
    if ("alpha_sigma2" %in% names(params)) {
      alpha_sigma2 <- params[["alpha_sigma2"]]
    }
    else {
      v <- mean(apply(Y, 1, var))
      alpha_sigma2 <- 2 + v / 4
      message(paste("alpha_sigma2 was not provided so I am using", alpha_sigma2))
    }
    if ("beta_sigma2" %in% names(params)) {
      beta_sigma2 <- params[["beta_sigma2"]]
    }
    else {
      if (!exists("v")) v <- mean(apply(Y, 1, var))
      beta_sigma2 <- (alpha_sigma2 - 1) * v
      message(paste("beta_sigma2 was not provided so I am using", beta_sigma2))
    }
    if (alpha_sigma2 <= 0) {
      stop("alpha_sigma2 is not positive; specify manually in params")
    }
    else if (beta_sigma2 <= 0) {
      stop("beta_sigma2 is not positive; specify manually in params")
    }
  }
  else if ("sigma2" %in% names(params)) {
    sigma2 <- params[["sigma2"]]
    if (!is.numeric(sigma2) || sigma2 <= 0) {
      stop("inavlid value for sigma2; must be numeric greater than zero")
    }
  }
  else {
    v <- mean(apply(Y, 1, var))
    sigma2 <- v
    message(paste("sigma2 was not provided so I am using", sigma2))
  }

  # Process Model; creates G_0 (mu_G) and Sigma_G_inv
  Sigma_G_inv <- matrix()
  if (proc_model == "RW") {
    G_0 <- diag(p)
  }
  else if (proc_model == "AR") {
    if ("mu_G" %in% names(params)) {
      if (is.diagonal.matrix(params$mu_G) &&
          all(dim(params$mu_G) != c(p, p))) {
        G_0 <- params[["mu_G"]]
      }
      else {
        stop("mu_G must be a diagonal p by p matrix")
      }

    }
    else {
      message("mu_G was not provided, so I am using an identity matrix")
      G_0 <- diag(p)
    }

    if ("Sigma_G_inv" %in% names(params)) {
      if (is.positive.definite(params$Sigma_G_inv) &&
          is.symmetric.matrix(params$Sigma_G_inv)) {
        Sigma_G_inv <- params[["Sigma_G_inv"]]
      } else {
        stop("Sigma_G_inv must be symmetric positive definite matrix")
      }
    } else {
      message("Sigma_G_inv was not provided, so I am using an identity matrix")
      Sigma_G_inv <- diag(p)
    }

  }
  else if (proc_model == "Full") {
    if ("mu_G" %in% names(params)) {
      if (is.matrix(params$mu_G) && all(dim(params$mu_G) != c(p, p))) {
        G_0 <- params$mu_G
      }
      else {
        stop("mu_G must be a p by p matrix")
      }
    }
    else {
      message("mu_G was not provided, so I am using an identity matrix")
      G_0 <- diag(p)
    }

    if ("Sigma_G_inv" %in% names(params)) {
      if (is.positive.definite(params$Sigma_G_inv) &&
          is.symmetric.matrix(params$Sigma_G_inv) &&
          nrow(params$Sigma_G_inv) == p^2) {
        Sigma_G_inv <- params$Sigma_G_inv
      }
      else {
        stop("Sigma_G_inv must be a p^2 by p^2 symmetric positive definite matrix")
      }
    }
    else {
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
      if (!is.numeric(alpha_lambda) || alpha_lambda <= 0) {
        stop("alpha_lambda must be numeric greater than 0")
      }
    }
    else {
      message("alpha_lambda was not provided so I am using 2.25")
      alpha_lambda <- 2.25
    }

    if ("beta_lambda" %in% names(params)) {
      beta_lambda <- params[["beta_lambda"]]
      if (!is.numeric(beta_lambda) || beta_lambda <= 0) {
        stop("beta_lambda must be numeric greater than 0")
      }
    }
    else {
      message("beta_lambda was not provided so I am using 0.0625")
      beta_lambda <- 0.0625
    }

    # Group scalar params into vector
    scalar_params <- c(alpha_lambda, beta_lambda, alpha_sigma2, beta_sigma2, sigma2)

    # Run the model
    results <- dstm_discount(Y, F_, G_0, Sigma_G_inv, m_0, C_0,
                  scalar_params, proc_model, n_samples, verbose)
  }
  else if (proc_error == "IW") {
    # C_W
    if ("C_W" %in% names(params)) {
      C_W <- params[["C_W"]]
      if (!is.positive.definite(C_W)) {
        stop("C_W must be a square positive definite matrix")
      }
    }
    else {
      message("C_W was not provided so I am using an identity matrix")
      C_W <- diag(p)
    }

    if ("df_W" %in% names(params)) {
      df_W <- as.integer(params[["df_W"]])
      if (!is.numeric(params[["df_W"]] || params[["df_W"]] < p)) {
        stop("df_W must be numeric >= p")
      }
    }
    else {
      message("df_W was not provided so I am using p")
      df_W <- p
    }

    scalar_params <- c(df_W, alpha_sigma2, beta_sigma2, sigma2)
    results <- dstm_IW(Y, F_, G_0, Sigma_G_inv, m_0, C_0, C_W,
                       scalar_params, proc_model, n_samples, verbose)
  }
  else {
    stop("I don't know that type of process error")
  }

  # Process output
  class(results) <- c("dstm", "list")
  attr(results,  "obs_model") <- obs_model
  attr(results, "proc_model") <- proc_model
  attr(results, "proc_error") <- proc_error
  attr(results, "sample_sigma2") <- sample_sigma2

  return(results)
}
