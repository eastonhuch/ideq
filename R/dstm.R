#' Fits a dynamic spatio-temporal model with EOFs
#'
#' @param Y S by T matrix containing response variable at S spatial locations and T time points
#' @param model character string; options include `discount`, `sample_G`, `AR`, and `IDE`
#' @param n_samples integer; number of posterior samples to take
#' @param P integer; dimension of G in the state equation \eqn{\theta_{t+1} = G \theta_{t}}
#' @param verbose boolean; controls verbosity
#' @param sample_sigma2 whether boolean; to sample \eqn{\sigma^2}
#'
#' @export
#' @examples
#' # Duhh...nothing yet
dstm_eof <- function(Y, proc_model = "RW",
                     proc_error = "IW", P = 10L,
                     n_samples = 1L, sample_sigma2 = TRUE,
                     verbose = FALSE, params = NULL) {

  # Observation Model; creates F, m_0, C_0
  e <- eigen(cov(t(Y)))
  F_ <-e$vectors[, 1:P]
  
  # Set m_0
  m_0 <- NULL
  if ("m_0" %in% names(params)) {
    if(length(params[["m_0"]]) != P) stop("m_0 must have length P")
    m_0 <- params[["m_0"]]
  } else {
    message("No prior was provided for m_0 so I am using a vector of zeros")
    m_0 <- rep(0, P)
  }

  # Set C_0
  C_0 <- matrix()
  if ("C_0" %in% names(params)) {
    if (!is.matrix(params[["C_0"]]) || any(dim(params[["C_0"]]) != c(P, P))){
      stop("C_0 must be a P by P matrix")
    }
    C_0 <- params[["C_0"]]
  }
  else {
    message("No prior was provided for C_0 so I am using 1e-6I")
    C_0 <- diag(1e-3, P)
  }

  # Observation Error; creates alpha_sigma2, beta_sigma2, sigma2
  # NOTE: We could also draw this from a Wishart distribution...nah
  alpha_sigma2 <- beta_sigma2 <- sigma2 <- -1
  if (sample_sigma2) {
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
    
    if (alpha_sigma2 <= 0) {
      stop("alpha_sigma2 is not positive; specify manually in params")
    }
    else if (beta_sigma2 <= 0) {
      stop("beta_sigma2 is not positive; specify manually in params")
    }
  } else if ("sigma2" %in% names(params)) {
    sigma2 <- params[["sigma2"]]
    if (!is.numeric(sigma2) || sigma2 <= 0) {
      stop("inavlid value for sigma2; must be numeric greater than zero")
    }
  } else {
    v <- mean(apply(Y, 1, var))
    sigma2 <- v
    message(paste("sigma2 was not provided so I am using", sigma2))
  }

  # Process Model; creates G_0 (mu_G) and Sigma_G_inv
  Sigma_G_inv <- matrix()
  if (proc_model == "RW") {
    G_0 <- diag(P)
  } else if (proc_model == "AR") {
    if ("mu_G" %in% names(params)) {
      if (matrixcalc::is.diagonal.matrix(params$mu_G) &&
          all(dim(params$mu_G) != c(P, P))) {
        G_0 <- params[["mu_G"]]
      }
      else {
        stop("mu_G must be a diagonal P by P matrix")
      }

    }
    else {
      message("mu_G was not provided, so I am using 10I")
      G_0 <- 1e1 * diag(P)
    }

    if ("Sigma_G_inv" %in% names(params)) {
      if (matrixcalc::is.positive.definite(params$Sigma_G_inv) &&
          matrixcalc::is.symmetric.matrix(params$Sigma_G_inv)) {
        Sigma_G_inv <- params[["Sigma_G_inv"]]
      } else {
        stop("Sigma_G_inv must be symmetric positive definite matrix")
      }
    } else {
      message("Sigma_G_inv was not provided, so I am using 1e3I")
      Sigma_G_inv <- diag(P)
    }

  } else if (proc_model == "Full") {
    if ("mu_G" %in% names(params)) {
      if (is.matrix(params$mu_G) && all(dim(params$mu_G) != c(P, P))) {
        G_0 <- params$mu_G
      }
      else {
        stop("mu_G must be a P by P matrix")
      }
    }
    else {
      message("mu_G was not provided, so I am using an identity matrix")
      G_0 <- diag(P)
    }

    if ("Sigma_G_inv" %in% names(params)) {
      if (matrixcalc::is.positive.definite(params$Sigma_G_inv) &&
          matrixcalc::is.symmetric.matrix(params$Sigma_G_inv) &&
          nrow(params$Sigma_G_inv) == P^2) {
        Sigma_G_inv <- params$Sigma_G_inv
      }
      else {
        stop("Sigma_G_inv must be a P^2 by P^2 symmetric positive definite matrix")
      }
    }
    else {
      message("Sigma_G_inv was not provided, so I am using 1e5I")
      Sigma_G_inv <- 1e5*diag(P^2)
    }

  } else {
    stop("proc_model not supported")
  }

  # Process Error; creates all necessary params (e.g., alpha_lambda)
  alpha_lambda <- beta_lambda <- df_W <- NA
  C_W <- matrix(NA)
  if (proc_error == "discount") {
    # Set prior for lambda
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
  } else if (proc_error == "IW") {
    # C_W
    if ("C_W" %in% names(params)) {
      C_W <- params[["C_W"]]
      if (!matrixcalc::is.positive.definite(C_W)) {
        stop("C_W must be a square positive definite matrix")
      }
    }
    else {
      message("C_W was not provided so I am using an identity matrix")
      C_W <- diag(P)
    }

    if ("df_W" %in% names(params)) {
      df_W <- as.integer(params[["df_W"]])
      if (!is.numeric(params[["df_W"]] || params[["df_W"]] < P)) {
        stop("df_W must be numeric >= P")
      }
    }
    else {
      message("df_W was not provided so I am using P")
      df_W <- P
    }
  } else {
    stop("I don't know that type of process error")
  }
  
  # Group scalar params into vector
  scalar_params <- c(alpha_lambda=alpha_lambda, beta_lambda=beta_lambda,
                     alpha_sigma2=alpha_sigma2, beta_sigma2=beta_sigma2,
                     sigma2=sigma2, df_W=df_W)

  # Fit the model
  results <- eof(Y, F_, G_0, Sigma_G_inv, m_0, C_0, C_W,
                 scalar_params, proc_model, n_samples, verbose)
  
  # Process results
  if ("lambda" %in% names(results)) {
    results[["lambda"]] <- as.numeric(results[["lambda"]])
  }
  if ("sigma2" %in% names(results)) {
    results[["sigma2"]] <- as.numeric(results[["sigma2"]])
  }

  # Process output
  class(results) <- c("dstm_eof", "dstm", "list")
  attr(results, "proc_model") <- proc_model
  attr(results, "proc_error") <- proc_error
  attr(results, "sample_sigma2") <- sample_sigma2

  return(results)
}

#' Fits an integrodifference equation model (IDE)
#' @description Detailed description
#'
#' @param Y matrix.
#'          S by T matrix containing response variable at S spatial locations and T time points
#' @param locs matrix or array.
#'             Spatial locations of observed data. If the locations are the same at each time
#'             point, then `locs` is an S by 2 matrix. 
#'             If the locations change at different time points,
#'             then `locs` is an array with dimension (S, 2, T)
#' @param kernel_locs integer or matrix.
#'                    if integer, then kernel parameters are calculated on a 2-D grid
#'                    with dimension (`kernel_locs`, `kernel_locs`)
#' @param proc_error "discount" or "IW".
#'                   If "discount", then the process error is estimated using a discount factor.
#'                   If "IW", then the process error is estimated as a realization from
#'                   an inverse-Wishart distribution
#' @param J integer. Extent of Fourier approximation.
#'          The size of the state space is 2 * J^2 + 1
#' @param n_samples integer; number of posterior samples to draw
#' @param sample_sigma2 boolean; whether to sample \eqn{\sigma^2}
#' @param verbose boolean; verbose = TRUE prints sampling progress
#' @param params list; Used to specify hyperparameters
#' 
#' @details How to specify hyperparameters such as priors
#'
#' @export
#' @examples
#' # Duhh...nothing yet
dstm_ide <- function(Y, locs=NULL, kernel_locs=NULL, proc_error = "discount", J=4L,
                     n_samples = 1L, sample_sigma2 = TRUE,
                     verbose = FALSE, params = NULL) {

  # Error checking for J, L, locs
  if (is.null(locs)) {
    stop("locs must be specified for obs_model == IDE")
  }

  J <- as.integer(J)
  if (!(is.numeric(J)) || J<1 ) {
    stop("J must an integer > 0")
  }

  if ("L" %in% names(params)) {
    L <- params[["L"]]
    L <- as.numeric(L)
    if (! (is.numeric(L) || L>0) ) {
      stop("L must be numeric > 0")
    }
  } else {
    L <- 2
    message(paste("L was not provided so I am using L="), L,
            " and centering and scaling locs")
    locs <- center_all(locs, L)
  }

  # Observation Model: creates m_0, C_0
  # Set m_0
  P <- 2*J^2 + 1
  m_0 <- NULL
  if ("m_0" %in% names(params)) {
    if(length(params[["m_0"]]) != P) stop("m_0 must have length P")
    m_0 <- params[["m_0"]]
  } else {
    message("No prior was provided for m_0 so I am using a vector of zeros")
    m_0 <- rep(0, P)
  }

  # Set C_0
  C_0 <- matrix()
  if ("C_0" %in% names(params)) {
    if (!is.matrix(params[["C_0"]]) || any(dim(params[["C_0"]]) != c(P, P))){
      stop("C_0 must be a P by P matrix")
    }
    C_0 <- params[["C_0"]]
  } else {
    message("No prior was provided for C_0 so I am using I/9")
    C_0 <- diag(1/9, P)
  }
  
  # Observation Error; creates alpha_sigma2, beta_sigma2, sigma2
  # NOTE: We could also draw this from a Wishart distribution
  alpha_sigma2 <- beta_sigma2 <- sigma2 <- -1
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
  } else if ("sigma2" %in% names(params)) {
    sigma2 <- params[["sigma2"]]
    if (!is.numeric(sigma2) || sigma2 <= 0) {
      stop("inavlid value for sigma2; must be numeric greater than zero")
    }
  } else {
    v <- mean(apply(Y, 1, var))
    sigma2 <- v
    message(paste("sigma2 was not provided so I am using", sigma2))
  }

  # Process Model; creates kernel parameters
  # FIXME: Add ability to use different locs
  
  num_kernel_locs <- 1
  locs_dim <- ncol(locs)
  if ("mu_kernel_mean" %in% names(params)) {
    mu_kernel_mean <- params[["mu_kernel_mean"]]
    mu_kernel_mean <- as.vector(mu_kernel_mean)
    if (!is.numeric(mu_kernel_mean) || length(mu_kernel_mean)!=locs_dim) {
      stop("mu_kernel_mean must be numeric with length==ncol(locs)")
    }
  } else {
    mu_kernel_mean <- rep(0, locs_dim)
    message("mu_kernel_mean was not provided so I am using a vector of zeroes")
  }
  
  if ("mu_kernel_var" %in% names(params)) {
    mu_kernel_var <- params[["mu_kernel_var"]]
    if (!matrixcalc::is.positive.definite(mu_kernel_var) ||
        any(dim(mu_kernel_var) != locs_dim)) {
      stop("mu_kernel_var must be positive definite with dimensions == ncol(locs)")
    }
  } else {
    mu_kernel_var <- diag(1/9, locs_dim)
    message("mu_kernel_var was not provided so I am using I/9")
  }
  
  if ("proposal_factor_mu" %in% names(params)) {
    proposal_factor_mu <- params[["proposal_factor_mu"]]
    if (!is.numeric(proposal_factor_mu) || proposal_factor_mu <= 0) {
      stop("proposal_factor_mu must be numeric > 0")
    }
  } else {
    proposal_factor_mu <- 1
    message(paste("proposal_factor_mu was not provided so I am using", 
            proposal_factor_mu))
  }
  
  if ("Sigma_kernel_df" %in% names(params)) {
    Sigma_kernel_df <- params[["Sigma_kernel_df"]]
    if (!is.numeric(Sigma_kernel_df) || Sigma_kernel_df <= 0) {
      stop("Sigma_kernel_df must be numeric > 0")
    }
  } else {
    Sigma_kernel_df <- 10
    message(paste("Sigma_kernel_df was not provided so I am using", 
            Sigma_kernel_df))
  }
  
  if ("Sigma_kernel_scale" %in% names(params)) {
    Sigma_kernel_scale <- params[["Sigma_kernel_scale"]]
    if (!matrixcalc::is.positive.definite(Sigma_kernel_scale) ||
        any(dim(Sigma) != locs_dim)) {
      stop("Sigma_kernel_scale must be positive definite with dimensions == ncol(locs)")
    }
  } else {
    Sigma_kernel_scale <- diag(1/Sigma_kernel_df, locs_dim)
    message("Sigma_kernel_scale was not provided so I am using I/Sigma_kernel_df")
  }
  
  if ("proposal_factor_Sigma" %in% names(params)) {
    proposal_factor_Sigma <- params[["proposal_factor_Sigma"]]
    if (!is.numeric(proposal_factor_Sigma) || proposal_factor_Sigma <= 0) {
      stop("proposal_factor_Sigma must be numeric > 0")
    }
  } else {
    proposal_factor_Sigma <- 1
    message(paste("proposal_factor_Sigma was not provided so I am using", 
            proposal_factor_Sigma))
  }
  
  # Error checking for kernel_locs
  # Also, adjustment to above quantities if using spatially varying kernel params
  SV <- FALSE
  K  <- array(0, dim=c(1, 1, 1))
  if (is.numeric(kernel_locs)) {
    SV <- TRUE
    # prepare kernel_locs
    if (length(kernel_locs) > 1) {
      kernel_locs <- center_all(kernel_locs, L)
    } else {
      kernel_locs <- gen_grid(as.integer(kernel_locs), L)
    }
    
    # Modify kernel parameters
    num_kernel_locs <- nrow(kernel_locs)
    mu_kernel_mean <- rep(1, num_kernel_locs) %x% mu_kernel_mean
    R <- exp(-as.matrix(dist(kernel_locs))) # Correlation matrix
    mu_kernel_var <- R %x% mu_kernel_var
    
    # Create K matrix
    K <- pdist::pdist(kernel_locs, locs)
    K <- as.matrix(K)
    K <- apply(K, 2, function(x) x / sum(x))
    K <- array(t(K), dim=c(dim(K), 1))
    
  } else if (is.null(kernel_locs)) {
    SV <- FALSE
  } else {
    stop("kernel_locs must be numeric or NULL")
  }
  
  if (is.null(num_kernel_locs)) num_kernel_locs <- 1
  Sigma_kernel_scale <- array(1, dim=c(1, 1, num_kernel_locs)) %x%
                        Sigma_kernel_scale
  
  # Process Error; creates alpha_lambda, beta_lambda, etc.
  alpha_lambda <- beta_lambda <- df_W <- NA
  C_W <- matrix(NA)
  if (proc_error == "discount") {
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
  } else if (proc_error == "IW") {
    # C_W
    if ("C_W" %in% names(params)) {
      C_W <- params[["C_W"]]
      if (!matrixcalc::is.positive.definite(C_W)) {
        stop("C_W must be a square positive definite matrix")
      }
    }
    else {
      message("C_W was not provided so I am using an identity matrix")
      C_W <- diag(P)
    }

    if ("df_W" %in% names(params)) {
      df_W <- as.integer(params[["df_W"]])
      if (!is.numeric(params[["df_W"]] || params[["df_W"]] < P)) {
        stop("df_W must be numeric >= P")
      }
    }
    else {
      message("df_W was not provided so I am using P")
      df_W <- P
    } 
  } else {
    stop("proc_error not recognized; must be `discount` or `IW`")
  }
  
  scalar_params <- c(alpha_sigma2=alpha_sigma2, beta_sigma2=beta_sigma2,
                     sigma2=sigma2, J=J, L=L, df_W=df_W,
                     alpha_lambda=alpha_lambda, beta_lambda=beta_lambda, 
                     proposal_factor_mu=proposal_factor_mu,
                     proposal_factor_Sigma=proposal_factor_Sigma,
                     Sigma_kernel_df=Sigma_kernel_df, SV=SV)
  
  results <- ide(Y, locs, m_0, C_0, mu_kernel_mean,
                 mu_kernel_var, K, Sigma_kernel_scale, C_W,
                 scalar_params, n_samples, verbose)
  
  # Process output
  class(results) <- c("dstm_ide" , "dstm", "list")
  attr(results, "proc_error") <- proc_error
  attr(results, "sample_sigma2") <- sample_sigma2

  return(results)
}
