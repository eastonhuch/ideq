#' Fits a dynamic spatio-temporal model (DSTM)
#'
#' @param Y S by T matrix containing response variable at S spatial locations and T time points
#' @param model character string; options include `discount`, `sample_G`, `AR`, and `IDE`
#' @param n_samples integer; number of posterior samples to take
#' @param p integer; dimension of G in the state equation \eqn{\theta_{t+1} = G \theta_{t}}
#' @param verbose boolean; controls verbosity
#' @param sample_sigma2 whether boolean; to sample \eqn{\sigma^2}
#'
#' @export
#' @examples
#' # Duhh...nothing yet
dstm <- function(Y, locs=NULL, obs_model = "EOF", proc_model = "RW",
                 proc_error = "discount", p = 10L,
                 n_samples = 1L, sample_sigma2 = TRUE,
                 verbose = FALSE, params = NULL) {
  results <- "No output...whoops"

  # Observation Model; creates F, m_0, C_0
  message("Observation model")
  IDE = obs_model == "IDE"
  if (IDE) {
    # Error checking for J, L, locs
    if (is.null(locs)) {
      stop("locs must be specified for obs_model == IDE")
    }

    if ("J" %in% names(params)) {
      J <- params[["J"]]
      J <- as.integer(J)
      if (! (is.numeric(J) || J>0) ) {
        stop("J must an integer > 0")
      }
    }
    else {
      J <- 5
      message(paste("J was not provided, so I am using "), J)
    }

    if ("L" %in% names(params)) {
      L <- params[["L"]]
      L <- as.numeric(L)
      if (! (is.numeric(L) || L>0) ) {
        stop("L must be numeric > 0")
      }
    }
    else {
      L <- 10
      message(paste("L was not provided so I am using L="), L,
              " and centering and scaling locs")
      center_col <- function(x, L) {
        x_range <- diff(range(x))
        x <- L/x_range * x
        x <- x + (L/2 - max(x))
        x
      }
      for ( i in seq(ncol(locs)) ) {
        locs[, i] <- center_col(locs[, i], L)
      }
    }

    # Set m_0
    p <- 2*J^2 + 1
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
      message("No prior was provided for C_0 so I am using 100I")
      C_0 <- 100*diag(p)
    }
  }
  else if (obs_model == "EOF") {
    e <- eigen(cov(t(Y)))
    F_ <-e$vectors[, 1:p]
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
      message("No prior was provided for C_0 so I am using diag(eigenvalues)")
      C_0 <- diag(e$values[1:p])
      #message("No prior was provided for C_0 so I am using I")
      #C_0 <- diag(p)
    }

    }
  else {
    stop("obs_model is invalid")
  }

  # Observation Error; creates alpha_sigma2, beta_sigma2, sigma2
  # NOTE: We could also draw this from a Wishart distribution...nah
  message("Observation error")
  alpha_sigma2 <- beta_sigma2 <- sigma2 <- NA
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
  message("Process model")
  Sigma_G_inv <- matrix()
  if (IDE) {
    # Error checking for m_kernel, C_kernel: the priors for the kernel parameters
    locs_dim <- ncol(locs)
    if ("m_kernel" %in% names(params)) {
      m_kernel <- params[["m_kernel"]]
      m_kernel <- as.vector(m_kernel)
      if (!is.numeric(m_kernel) || length(m_kernel)!=locs_dim) {
        stop("m_kernel must be numeric with length==ncol(locs)")
      }
    }
    else {
      m_kernel <- rep(0, locs_dim)
      message("m_kernel was not provided so I am using a vector of zeroes")
    }

    if ("C_kernel" %in% names(params)) {
      C_kernel <- params[["C_kernel"]]
      if (!is.positive.definite(C_kernel) ||
          any(dim(Sigma) != locs_dim)) {
        stop("C_kernel must be positive definite with dimensions == ncol(locs)")
      }
    }
    else {
      C_kernel <- diag(locs_dim)
      message("C_kernel was not provided so I am using I")
    }
  }
  else if (proc_model == "RW") {
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
      message("mu_G was not provided, so I am using I")
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
      message("Sigma_G_inv was not provided, so I am using 100I")
      Sigma_G_inv <- 100*diag(p)
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
      message("Sigma_G_inv was not provided, so I am using 100I")
      Sigma_G_inv <- 100*diag(p^2)
    }

  }
  else {
    stop("proc_model not supported")
  }

  # Process Error; creates all necessary params (e.g., alpha_lambda)
  message("Process error")
  if (IDE) {
    # FIXME: Add code to to specify process error in IDE model
    if ("alpha_lambda" %in% names(params)) {
      alpha_lambda <- params[["alpha_lambda"]]
    }
    else {
      alpha_lambda <- 4
      message(paste("alpha_lambda was not provided so I am using", alpha_lambda))
    }
    if ("beta_lambda" %in% names(params)) {
      beta_lambda <- params[["beta_lambda"]]
    }
    else {
      beta_lambda <- 3
      message(paste("beta_lambda was not provided so I am using", beta_lambda))
    }

    scalar_params <- c(alpha_sigma2=alpha_sigma2, beta_sigma2=beta_sigma2,
                       sigma2=sigma2, J=J, L=L,
                       alpha_lambda=alpha_lambda, beta_lambda=beta_lambda)
    results <- dstm_IDE(Y, locs, m_0, C_0, m_kernel, C_kernel,
                        scalar_params, n_samples, verbose)
  }
  else if (proc_error == "discount") {

    # Set prior for lambda
    alpha_lambda <- beta_lambda <- NULL
    if ("alpha_lambda" %in% names(params)) {
      alpha_lambda <- params[["alpha_lambda"]]
      if (!is.numeric(alpha_lambda) || alpha_lambda <= 0) {
        stop("alpha_lambda must be numeric greater than 0")
      }
    }
    else {
      alpha_lambda <- 3
      message(paste("alpha_lambda was not provided so I am using", alpha_lambda))
    }

    if ("beta_lambda" %in% names(params)) {
      beta_lambda <- params[["beta_lambda"]]
      if (!is.numeric(beta_lambda) || beta_lambda <= 0) {
        stop("beta_lambda must be numeric greater than 0")
      }
    }
    else {
      beta_lambda <- 4
      message(paste("beta_lambda was not provided so I am using", beta_lambda))
    }

    # Group scalar params into vector
    scalar_params <- c(alpha_lambda=alpha_lambda, beta_lambda=beta_lambda,
                       alpha_sigma2=alpha_sigma2, beta_sigma2=beta_sigma2,
                       sigma2=sigma2)

    # Run the model
    results <- dstm_discount(Y, F_, G_0, Sigma_G_inv, m_0, C_0,
                  scalar_params, proc_model, n_samples, verbose)
    results[["lambda"]] <- as.numeric(results[["lambda"]])
    if ("sigma2" %in% names(results)) {
      results[["sigma2"]] <- as.numeric(results[["sigma2"]])
    }
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

    scalar_params <- c(df_W=df_W, sigma2=sigma2,
                       alpha_sigma2=alpha_sigma2, beta_sigma2=beta_sigma2)
    results <- dstm_IW(Y, F_, G_0, Sigma_G_inv, m_0, C_0, C_W,
                       scalar_params, proc_model, n_samples, verbose)
    if ("sigma2" %in% names(results)) {
      results[["sigma2"]] <- as.numeric(results[["sigma2"]])
    }
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
