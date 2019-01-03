#' Dynamic spatio-temporal model with EOFs
#' 
#' @description 
#' Fits a dynamic spatio-temporal model using empirical orthogonal functions
#' (EOFs).
#' The model does not require the spatial locations because the process model
#' is based on the principal components of the data matrix.
#' Three broad model types are supported:
#' 
#' 1. RW: A random walk model for which the process matrix is the identity.
#' 
#' 2. AR: An auto-regressive model for which the process matrix is diagonal
#'    and its elements are estimated.
#'    
#' 3. Dense: A model in which the process matrix is a dense, estimated
#'    matrix.
#'    
#' For each broad model type,
#' users can specify a variety of options including
#' the size of the state space,
#' the form of the process error,
#' and whether to sample the observation error.
#' Users can specify prior distributions for all sampled quantities using the
#' `params` argument.
#' 
#' @param Y 
#' (numeric matrix) S by T data matrix containing the response variable
#' at S spatial locations and T time points.
#' The t-th column (NOT row) corresponds to the t-th observation vector.
#' @param proc_model 
#' (character string) Process model: one of "RW" (identity process matrix),
#' "AR" (diagonal process matrix), or "Dense" (dense process matrix).
#' @param P 
#' (integer) Number of EOFs.
#' @param proc_error 
#' (character string) Process error: 
#' "IW" (inverse-Wishart) or "Discount" (discount factor).
#' @param sample_sigma2 
#' (logical) whether to sample the variance of the iid observation error.
#' @param verbose 
#' (logical) Whether to print additional information;
#' e.g., iteration in sampling algorithm.
#' @param params 
#' (list) List of hyperparameter values; see details.
#' 
#' @details 
#' This section explains how to specify custom hyperparameters using the `params` argument.
#' For each distribution referenced below,
#' we use the scale parameterization found on the distribution's Wikipedia page.
#' You may specify the following as named elements of the `params` list:
#' 
#' m_0: (numeric vector) The prior mean of the state vector at time zero
#'  (\eqn{\theta_0})
#'
#' C_0: (numeric matrix) The prior variance-covariance matrix of the state
#' vector at time zero (\eqn{\theta_0})
#'
#' alpha_sigma2, beta_sigma2: (numeric scalars) The inverse-Gamma parameters 
#' (scale parameterization) of the prior distribution on \eqn{\sigma^2}
#' 
#' sigma2: (numeric scalar) The value to use for \eqn{\sigma^2}
#' if sample_sigma2 = FALSE
#' 
#' mu_G: (numeric matrix) The prior mean for the process matrix G.
#' If proc_model = "AR", then mu_G must also be diagonal.
#' If proc_model = "Dense", then mu_G has no contraints.
#' 
#' Sigma_G: (numeric matrix) The prior variance-covariance matrix for the
#' process matrix. If proc_model = "AR", then Sigma_G should be P by 
#' P and is the variance-covariance matrix for diag(G).
#' If proc_model = "Dense", then Sigma_G should be P^2 by P^2 and is the 
#' variance-covariance matrix for vec(G)
#' 
#' alpha_lambda, beta_lambda: (numeric scalars) The inverse-Gamma parameters 
#' (scale parameterization) of the prior distribution on 
#' \eqn{\lambda = \delta / (1 - \delta)}
#' 
#' C_W: (numeric matrix) The scale matrix for the inverse-Wishart prior
#' distribution on W, the variance-covariance matrix of the process error.
#' 
#' df_W: (numeric scalar) The degees of freedom for the inverse-Wishart prior
#' distribution on W, the variance-covariance matrix of the process error.
#' 
#' @examples
#' # Create example data
#' num_time_points <- 5
#' num_spatial_locations <- 100
#' z <- rnorm(num_time_points * num_spatial_locations)
#' Y <- matrix(z, nrow=num_spatial_locations, ncol=num_time_points)
#' 
#' # Illustrate methods
#' rw_model <- dstm_eof(Y, proc_model="RW") # Random walk model
#' summary(rw_model) # print(rw_model) is equivalent
#' predictions <- predict(rw_model) 
#' 
#' # Other model types
#' dstm_eof(Y, proc_model="AR") # Diagonal process matrix
#' dstm_eof(Y, proc_model="Dense") # Dense process matrix
#' dstm_eof(Y, proc_error="Discount") # Discount factor process error
#' 
#' # Specify hyperparameters
#' dstm_eof(Y, sample_sigma2=FALSE, params=list(sigma2=1)) # Fix sigma2
#' dstm_eof(Y, params=list(alpha_lambda=10, beta_lambda=11)) # Prior for lambda
#' dstm_eof(Y, P=10, params=list(m_0=rep(1, 10) , C_0=diag(0.01, 10))) # Prior for theta_0
#' dstm_eof(Y, params=list(C_W=diag(10), df_W=100)) # Prior for W
#' @export
dstm_eof <- function(Y, proc_model = "Dense", P = 10L, proc_error = "IW",
                     n_samples = 1L, sample_sigma2 = TRUE, verbose = FALSE,
                     params = NULL) {

  # Observation Model; creates F, m_0, C_0
  e <- eigen(cov(t(Y)))
  F_ <-e$vectors[, 1:P]
  
  # Set m_0
  m_0 <- NULL
  if ("m_0" %in% names(params)) {
    if(length(params[["m_0"]]) != P) stop("m_0 must have length P")
    m_0 <- params[["m_0"]]
  } else {
    message("m_0 was not provided so I am using a vector of zeros")
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
    k <- 1e-6
    message(paste("C_0 was not provided so I am using", k, "I"))
    C_0 <- diag(k, P)
  }

  # Observation Error; creates alpha_sigma2, beta_sigma2, sigma2
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
    
    if ( !is.numeric(alpha_sigma2) || alpha_sigma2 <= 0 ) {
      stop("alpha_sigma2 is not positive numeric; specify manually in params")
    }
    if ( !is.numeric(beta_sigma2) || beta_sigma2 <= 0 ) {
      stop("beta_sigma2 is not positive numeric; specify manually in params")
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

  # Process Model; creates G_0 (mu_G) and Sigma_G
  Sigma_G <- Sigma_G_inv <- matrix()
  if (proc_model == "RW") {
    G_0 <- diag(P)
  } else if (proc_model == "AR") {
    if ("mu_G" %in% names(params)) {
      if (matrixcalc::is.diagonal.matrix(params$mu_G) &&
          all(dim(params$mu_G) != P)) {
        G_0 <- params[["mu_G"]]
      }
      else {
        stop("mu_G must be a diagonal P by P matrix")
      }

    }
    else {
      message("mu_G was not provided, so I am using I")
      G_0 <- diag(P)
    }

    if ("Sigma_G" %in% names(params)) {
      if (matrixcalc::is.positive.definite(params$Sigma_G) &&
          matrixcalc::is.symmetric.matrix(params$Sigma_G) &&
          all(dim(params$Sigma_G) == P)) {
        Sigma_G <- params[["Sigma_G"]]
      } else {
        stop("Sigma_G must be P by P symmetric positive definite matrix")
      }
    } else {
      k <- 1e-3
      message(paste("Sigma_G was not provided, so I am using", k, "* I"))
      Sigma_G <- diag(k, P)
    }
    Sigma_G_inv <- chol2inv(chol(Sigma_G))

  } else if (proc_model == "Dense") {
    if ("mu_G" %in% names(params)) {
      if (is.matrix(params$mu_G) && all(dim(params$mu_G) == P)) {
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

    if ("Sigma_G" %in% names(params)) {
      if (matrixcalc::is.positive.definite(params$Sigma_G) &&
          matrixcalc::is.symmetric.matrix(params$Sigma_G) &&
          nrow(params$Sigma_G) == P^2) {
        Sigma_G <- params$Sigma_G
      }
      else {
        stop("Sigma_G must be a P^2 by P^2 symmetric positive definite matrix")
      }
    }
    else {
      k <- 1e-5
      message(paste("Sigma_G was not provided, so I am using", k, "* I"))
      Sigma_G <- diag(k, P^2)
    }
    Sigma_G_inv <- chol2inv(chol(Sigma_G))
    
  } else {
    stop("proc_model not supported")
  }

  # Process Error; creates all necessary params (e.g., alpha_lambda)
  alpha_lambda <- beta_lambda <- df_W <- NA
  C_W <- matrix(NA)
  if (proc_error == "Discount") {
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
      df_W <- as.numeric(params[["df_W"]])
      if (!is.numeric(P) || df_W < P) {
        stop("df_W must be numeric >= P")
      }
    }
    else {
      k <- 1e3
      message(paste("df_W was not provided so I am using", k, "* P"))
      df_W <- k * P
    }
  } else {
    stop("proc_error not recognized")
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

#' Integrodifference equation (IDE) model
#' 
#' @description
#' dstm_ide fits a type of dynamic spatio-temporal model called
#' an integrodifference equation (IDE) model.
#' It estimates a redistribution kernel---a
#' probability distribution controlling dispersion across time and space.
#' Currently, only normal redistribution kernels are supported.
#' 
#' The process model is decomposed with an orthonormal basis function expansion
#' (a Fourier series).
#' It can then be estimated as a special case of a dynamic linear model (DLM),
#' using the forward filtering backward sampling algorithm to estimate the
#' state vector (the Fourier coefficients).
#' The kernel parameters are estimated with a random walk Metropolis-Hastings
#' update.
#' The other parameters are estimated with conditionally conjugate updates.
#' 
#' @param Y 
#' (numeric matrix) S by T data matrix containing response variable at S spatial
#' locations and T time points.
#' The t-th column (NOT row) corresponds to the t-th observation vector.
#' @param locs
#' (numeric matrix)
#' S by 2 matrix containing the spatial locations of the observed data.
#' The rows of `locs` correspond with the rows of `Y`.
#' @param knot_locs 
#' (integer or numeric matrix) Knot locations.
#' The kernel parameters are estimated at these locations and then mapped to the
#' spatial locations of the observed data via process convolution.
#' If an integer is provided, then the knots are located on an equally spaced
#' grid with dimension (`knot_locs`, `knot_locs`).
#' If a matrix is provided, 
#' then each row of the matrix corresponds to a knot location.
#' @param proc_error 
#' (character string) Process error:
#' "IW" (inverse-Wishart) or "Discount" (discount factor).
#' @param J 
#' (integer) Extent of the Fourier approximation.
#' The size of the state space is 2 * J^2 + 1.
#' @param n_samples 
#' (integer) Number of posterior samples to draw.
#' @param sample_sigma2 
#' (logical) Whether to sample the variance of the iid observation error.
#' @param verbose 
#' (logical) Whether to print additional information;
#' e.g., iteration in sampling algorithm.
#' @param params 
#' (list) List of hyperparameter values; see details.
#' 
#' @details
#' This section explains how to specify custom hyperparameters using the `params` argument.
#' For each distribution referenced below,
#' we use the scale parameterization found on the distribution's Wikipedia page.
#' You may specify the following as named elements of the `params` list:
#' 
#' m_0: (numeric vector) The prior mean of the state vector
#' (Fourier basis coefficients) at time zero (\eqn{\theta_0}).
#'
#' C_0: (numeric matrix) The prior variance-covariance matrix of the state
#'  vector (Fourier basis coefficients) at time zero (\eqn{\theta_0}).
#'
#' alpha_sigma2, beta_sigma2: (numeric scalars) The inverse-Gamma parameters
#' (scale parameterization) of the prior distribution on \eqn{\sigma^2}.
#' 
#' sigma2: (numeric scalar) The value to use for \eqn{\sigma^2} 
#' if sample_sigma2 = FALSE.
#' 
#' alpha_lambda, beta_lambda: (numeric scalars) The inverse-Gamma parameters 
#' (scale parameterization) of the prior distribution on 
#' \eqn{\lambda = \delta / (1 - \delta)}.
#' 
#' C_W: (numeric matrix) The scale matrix for the inverse-Wishart prior
#' distribution on W, the variance-covariance matrix of the process error.
#' 
#' df_W: (numeric scalar) The degees of freedom for the inverse-Wishart prior
#' distribution on W, the variance-covariance matrix of the process error.
#' 
#' L: (numeric scalar) The range of the spatial locations after rescaling.
#' For spatially varying kernels, the value of L controls the degree of 
#' smoothing. As L increases, the degree of smoothing decreases.
#' 
#' mu_kernel_mean: (numeric vector) The mean of the normal prior distribution
#' on mu_kernel, the mean of the redistribution kernel.
#' In the spatially varying case, the prior distribution for mu_kernel
#' is assumed to be the same at every knot location.
#' 
#' mu_kernel_var: (numeric matrix) The variance of the normal prior distribution
#' on mu_kernel, the mean of the redistribution kernel.
#' 
#' Sigma_kernel_scale: (numeric matrix) The scale matrix for the 
#' inverse-Wishart prior distribution on Sigma_kernel,
#' the variance-covariance matrix of the redistribution kernel.
#' 
#' Sigma_kernel_df: (numeric scalar) The degrees of freedom for the 
#' inverse-Wishart prior distribution on Sigma_kernel,
#' the variance-covariance matrix of the redistribution kernel.
#' 
#' proposal_factor_mu: (numeric scalar) Controls the variance of the proposal distribution for
#' mu. The proposals have a variance of proposal_factor_mu^2 * mu_kernel_var.
#' proposal_factor_mu must generally be set lower for spatially varying models.
#' 
#' proposal_factor_Sigma: (numeric scalar) Controls the variance of the proposal distribution
#' for Sigma. As is the case with proposal_factor_mu, a higher value
#' corresponds to a higher variance.
#' The degrees of freedom for the proposal distribution for Sigma is 
#' ncol(locs) + Sigma_kernel_df / proposal_factor_Sigma.
#' proposal_factor_Sigma must generally be set lower for spatially varying 
#' models.
#' 
#'
#' @examples
#' # Create example data
#' num_time_points <- 5
#' spatial_locations <- expand.grid(seq(10), seq(10))
#' num_spatial_locations <- nrow(spatial_locations)
#' z <- rnorm(num_time_points * num_spatial_locations)
#' Y <- matrix(z, nrow=num_spatial_locations, ncol=num_time_points)
#' 
#' # Basic IDE model with one kernel
#' dstm_ide(Y, spatial_locations)
#' 
#' # IDE model with spatially varying kernel
#' dstm_ide(Y, spatial_locations, knot_locs=4)
#' 
#' # Discount factor method for estimating process error variance
#' dstm_ide(Y, spatial_locations, proc_error="Discount")
#' 
#' # Fix sigma2
#' dstm_ide(Y, spatial_locations, sample_sigma2=FALSE, params=list(sigma2=1))
#' 
#' # Set prior on sigma2
#' dstm_ide(Y, spatial_locations, 
#'          params=list(alpha_sigma2=10, beta_sigma2=11))
#' 
#' # Rescale spatial locations to have range of 10
#' dstm_ide(Y, spatial_locations, params=list(L=10)) 
#' 
#' # Set prior on kernel mean
#' dstm_ide(Y, spatial_locations, 
#'          params=list(mu_kernel_mean=c(0.2, 0.4),
#'                      mu_kernel_var=diag(2))) 
#' 
#' # Set prior on kernel variance-covariance matrix
#' dstm_ide(Y, spatial_locations, 
#'          params=list(Sigma_kernel_scale=diag(2), Sigma_kernel_df=100))
#' 
#' # Set prior on state vector (Fourier basis coefficients)
#' dstm_ide(Y, spatial_locations, 
#'          params=list(m_0=rep(0.1, 33), C_0=diag(0.01, 33))) 
#' 
#' # Set prior on process error
#' dstm_ide(Y, spatial_locations, 
#'          params=list(C_W=diag(33), df_W=100))
#' 
#' # Set proposal scaling factors
#' dstm_ide(Y, spatial_locations, 
#'          params=list(proposal_factor_mu=2,
#'                      proposal_factor_Sigma=3))
#' @export
dstm_ide <- function(Y, locs=NULL, knot_locs=NULL, proc_error = "IW", J=4L,
                     n_samples = 10L, sample_sigma2 = TRUE,
                     verbose = FALSE, params = NULL) {

  # Error checking for J, L, locs
  if (is.null(locs)) {
    stop("locs must be specified for IDE models")
  } else {
    if (class(locs) == "data.frame") locs <- as.matrix(locs)
    if (class(locs) != "matrix") {
      stop("locs must be data.frame or matrix")
    }
  }

  J <- as.integer(J)
  if ( is.na(J) || J<1 ) {
    stop("J must an integer > 0")
  }

  if ("L" %in% names(params)) {
    L <- params[["L"]]
    L <- as.numeric(L)
    if ( is.na(L) || L<=0 ) {
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
    message("m_0 was not provided so I am using a vector of zeros")
    m_0 <- rep(0, P)
  }

  # Set C_0
  C_0 <- matrix()
  if ("C_0" %in% names(params)) {
    C_0 <- params[["C_0"]]
    if (!is.matrix(C_0) || any(dim(C_0) != P)) {
      stop("C_0 must be a P by P matrix")
    }
  } else {
    message("C_0 was not provided so I am using I/9")
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
    if ( !is.numeric(alpha_sigma2) || alpha_sigma2 <= 0 ) {
      stop("alpha_sigma2 is not positive numeric; specify manually in params")
    }
    if ( !is.numeric(beta_sigma2) || beta_sigma2 <= 0 ) {
      stop("beta_sigma2 is not positive numeric; specify manually in params")
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
  
  n_knots <- 1
  locs_dim <- ncol(locs)
  if ("mu_kernel_mean" %in% names(params)) {
    mu_kernel_mean <- params[["mu_kernel_mean"]]
    mu_kernel_mean <- as.matrix(mu_kernel_mean)
    if (!is.numeric(mu_kernel_mean) || length(mu_kernel_mean)!=locs_dim) {
      stop("mu_kernel_mean must be numeric with length==ncol(locs)")
    }
  } else {
    mu_kernel_mean <- matrix(rep(0, locs_dim))
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
        any(dim(Sigma_kernel_scale) != locs_dim)) {
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
  
  # Error checking for knot_locs
  # Also, adjustment to above quantities if using spatially varying kernel params
  SV <- FALSE
  K  <- array(0, dim=c(1, 1, 1))
  if (is.numeric(knot_locs)) {
    SV <- TRUE
    # prepare knot_locs
    if (length(knot_locs) > 1) {
      knot_locs <- center_all(knot_locs, L)
    } else {
      knot_locs <- gen_grid(as.integer(knot_locs), L)
    }
    
    # Modify kernel parameters
    n_knots <- nrow(knot_locs)
    mu_kernel_mean <- rep(1, n_knots) %x% matrix(mu_kernel_mean, nrow=1)
    mu_kernel_var <- mu_kernel_var %x% diag(n_knots)
    
    # Create K matrix
    K <- pdist::pdist(knot_locs, locs)
    K <- as.matrix(K)
    K <- exp(-K)
    K <- apply(K, 2, function(x) x / sum(x))
    K <- t(K) # Makes K of dimension n_locs by n_knots
    K <- array(K, dim=c(dim(K), 1))
    
  } else if (is.null(knot_locs)) {
    SV <- FALSE
  } else {
    stop("knot_locs must be numeric or NULL")
  }
  
  if (is.null(n_knots)) n_knots <- 1
  Sigma_kernel_scale <- array(1, dim=c(1, 1, n_knots)) %x%
                        Sigma_kernel_scale
  
  # Process Error; creates alpha_lambda, beta_lambda, etc.
  alpha_lambda <- beta_lambda <- df_W <- NA
  C_W <- matrix(NA)
  if (proc_error == "Discount") {
    if ("alpha_lambda" %in% names(params)) {
      alpha_lambda <- params[["alpha_lambda"]]
      if ( !is.numeric(alpha_lambda) || alpha_lambda < 0 ) {
        stop("alpha_lambda must be numeric > 0")
      }
    }
    else {
      alpha_lambda <- 4
      message(paste("alpha_lambda was not provided so I am using", alpha_lambda))
    }
    
    if ("beta_lambda" %in% names(params)) {
      beta_lambda <- params[["beta_lambda"]]
      if ( !is.numeric(beta_lambda) || beta_lambda < 0 ) {
        stop("beta_lambda must be numeric > 0")
      }
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
        stop("C_W must be a positive definite matrix")
      }
    }
    else {
      message("C_W was not provided so I am using an identity matrix")
      C_W <- diag(P)
    }

    if ("df_W" %in% names(params)) {
      df_W <- params[["df_W"]]
      if ( !is.numeric(df_W) || df_W < P ) {
        stop("df_W must be numeric >= P")
      }
    }
    else {
      message("df_W was not provided so I am using P")
      df_W <- P
    } 
  } else {
    stop("proc_error not recognized; must be `Discount` or `IW`")
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
  
  # Process results
  if ("lambda" %in% names(results)) {
    results[["lambda"]] <- as.numeric(results[["lambda"]])
  }
  if ("sigma2" %in% names(results)) {
    results[["sigma2"]] <- as.numeric(results[["sigma2"]])
  }
  if (length(results[["Sigma_kernel"]]) > 1) {
    results[["Sigma_kernel"]] <- results[["Sigma_kernel"]][-1]
  }
  results[["mu_acceptance_rate"]] <- results[["mu_acceptances"]] / n_samples
  results[["Sigma_acceptance_rate"]] <- results[["Sigma_acceptances"]] / n_samples
  results[["mu_acceptances"]] <- NULL
  results[["Sigma_acceptances"]] <- NULL
  
  class(results) <- c("dstm_ide" , "dstm", "list")
  attr(results, "proc_model") <- "ide"
  attr(results, "proc_error") <- proc_error
  attr(results, "sample_sigma2") <- sample_sigma2

  return(results)
}
