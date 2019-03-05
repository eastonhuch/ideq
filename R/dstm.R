#' Dynamic spatio-temporal model with EOFs
#' @useDynLib ideq
#' @importFrom Rcpp sourceCpp
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
#' (integer) Number of EOFs or, in other words, the state space size.
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
#'  (\eqn{\theta_0}).
#'
#' C_0: (numeric matrix) The prior variance-covariance matrix of the state
#' vector at time zero (\eqn{\theta_0}).
#'
#' alpha_sigma2, beta_sigma2: (numeric scalars) The inverse-Gamma parameters 
#' (scale parameterization) of the prior distribution on the observation error 
#' (\eqn{\sigma^2}).
#' 
#' sigma2: (numeric scalar) The value to use for the observation error 
#' (\eqn{\sigma^2}) if sample_sigma2 = FALSE.
#' 
#' mu_G: (numeric matrix) The prior mean for the process matrix G.
#' If proc_model = "AR", then mu_G must be a diagonal matrix.
#' If proc_model = "Dense", then mu_G has no constraints.
#' 
#' Sigma_G: (numeric matrix) The prior variance-covariance matrix for the
#' process matrix. If proc_model = "AR", then Sigma_G should be P by 
#' P and is the variance-covariance matrix for diag(G).
#' If proc_model = "Dense", then Sigma_G should be P^2 by P^2 and is the 
#' variance-covariance matrix for vec(G).
#' 
#' alpha_lambda, beta_lambda: (numeric scalars) The inverse-Gamma parameters 
#' (scale parameterization) of the prior distribution on 
#' \eqn{\lambda = (1 - \delta) / \delta},
#' where \eqn{\delta} is the discount factor.
#' 
#' C_W: (numeric matrix) The scale matrix for the inverse-Wishart prior
#' distribution on the variance-covariance matrix of the process error (W).
#' 
#' df_W: (numeric scalar) The degees of freedom for the inverse-Wishart prior
#' distribution on the variance-covariance matrix of the process error (W).
#' 
#' @examples
#' # Load example data
#' data("standard_ide_data")
#' 
#' # Illustrate methods
#' rw_model <- dstm_eof(standard_ide_data, proc_model="RW", verbose=TRUE)
#' summary(rw_model) # print(rw_model) is equivalent
#' predictions <- predict(rw_model) 
#' 
#' # Other model types
#' dstm_eof(standard_ide_data, proc_model="AR") # Diagonal process matrix
#' dstm_eof(standard_ide_data, proc_model="Dense") # Dense process matrix
#' dstm_eof(standard_ide_data, proc_error="Discount") # Discount factor
#' 
#' # Specify hyperparameters
#' dstm_eof(standard_ide_data, sample_sigma2=FALSE, params=list(sigma2=0.01))
#' dstm_eof(standard_ide_data, proc_error="Discount", 
#'          params=list(alpha_lambda=201, beta_lambda=20))
#' dstm_eof(standard_ide_data, P=10, 
#'          params=list(m_0=rep(1, 10), C_0=diag(0.01, 10)))
#' dstm_eof(standard_ide_data, params=list(C_W=diag(10), df_W=100))
#' 
#' @export
dstm_eof <- function(Y, proc_model = "Dense", P = 10L, proc_error = "IW",
                     n_samples = 10L, sample_sigma2 = TRUE, verbose = FALSE,
                     params = NULL) {
  
  # Observation Model; creates F, m_0, C_0
  e <- eigen(cov(t(Y)))
  F_ <-e$vectors[, 1:P]
  
  # Process Model; creates mu_G and Sigma_G
  mu_G <- Sigma_G <- Sigma_G_inv <- matrix()
  
  # mu_G
  mu_G <- if ( is.null(params[["mu_G"]]) ) diag(P) else params[["mu_G"]]
  if (proc_model=="AR" && !matrixcalc::is.diagonal.matrix(mu_G))
    stop("mu_G must be diagonal for proc_model=\"AR\"")
  check.numeric.matrix(mu_G, P)
  
  # Sigma_G
  if (proc_model == "AR") {
    Sigma_G <- params[["Sigma_G"]] %else% diag(1e-2, P)
    check.dim(Sigma_G, P, "Sigma_G")
    check.cov.matrix(Sigma_G, P)
  } else if (proc_model == "Dense") {
    Sigma_G <- params[["Sigma_G"]] %else% diag(1e-2, P^2)
    check.dim(Sigma_G, P^2, "Sigma_G", "P^2")
    check.cov.matrix(Sigma_G, P^2)
  }
    
  # Sigma_G_inv
  Sigma_G_inv <- chol_inv(Sigma_G)
  
  # Process remaining params
  new_params <- process_common_params(params, proc_error, P, sample_sigma2)
  scalar_params <- c(alpha_lambda=new_params[["alpha_lambda"]], 
                     beta_lambda=new_params[["beta_lambda"]],
                     alpha_sigma2=new_params[["alpha_sigma2"]], 
                     beta_sigma2=new_params[["beta_sigma2"]],
                     sigma2=new_params[["sigma2"]], 
                     df_W=new_params[["df_W"]])

  # Fit the model
  results <- eof(Y, F_, mu_G, Sigma_G_inv, new_params[["m_0"]], 
                 new_params[["C_0"]], new_params[["C_W"]],
                 scalar_params, proc_model, n_samples, verbose)
  
  # Process results
  if ("lambda" %in% names(results))
    results[["lambda"]] <- as.numeric(results[["lambda"]])
  if ("sigma2" %in% names(results))
    results[["sigma2"]] <- as.numeric(results[["sigma2"]])

  # Process output
  class(results) <- c("dstm_eof", "dstm", "list")
  attr(results, "proc_model") <- proc_model
  attr(results, "proc_error") <- proc_error
  attr(results, "sample_sigma2") <- sample_sigma2

  return(results)
}

#' Integrodifference equation (IDE) model
#' @importFrom Rcpp sourceCpp
#' @description
#' dstm_ide fits a type of dynamic spatio-temporal model called
#' an integrodifference equation (IDE) model.
#' It estimates a redistribution kernel---a
#' probability distribution controlling diffusion across time and space.
#' Currently, only normal redistribution kernels are supported.
#' 
#' The process model is decomposed with an orthonormal basis function expansion
#' (a Fourier series).
#' It can then be estimated as a special case of a dynamic linear model (DLM),
#' using the forward filtering backward sampling algorithm to estimate the
#' state vector.
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
#' (integer or numeric matrix) Knot locations for the spatially varying IDE model.
#' The kernel parameters are estimated at these locations and then mapped to the
#' spatial locations of the observed data via process convolution.
#' If an integer is provided, then the knots are located on an equally spaced
#' grid with dimension (`knot_locs`, `knot_locs`).
#' If a matrix is provided, 
#' then each row of the matrix corresponds to a knot location.
#' If NULL, then the standard (spatially constant) IDE is fit.
#' @param proc_error 
#' (character string) Process error:
#' "IW" (inverse-Wishart) or "Discount" (discount factor).
#' @param J 
#' (integer) Extent of the Fourier approximation.
#' The size of the state space is (2*J + 1)^2.
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
#' m_0: (numeric vector) The prior mean of the state vector at time zero
#'  (\eqn{\theta_0}).
#'
#' C_0: (numeric matrix) The prior variance-covariance matrix of the state
#' vector at time zero (\eqn{\theta_0}).
#'
#' alpha_sigma2, beta_sigma2: (numeric scalars) The inverse-Gamma parameters 
#' (scale parameterization) of the prior distribution on the observation error 
#' (\eqn{\sigma^2}).
#' 
#' sigma2: (numeric scalar) The value to use for the observation error 
#' (\eqn{\sigma^2}) if sample_sigma2 = FALSE.
#' 
#' alpha_lambda, beta_lambda: (numeric scalars) The inverse-Gamma parameters 
#' (scale parameterization) of the prior distribution on 
#' \eqn{\lambda = (1 - \delta) / \delta},
#' where \eqn{\delta} is the discount factor.
#' 
#' C_W: (numeric matrix) The scale matrix for the inverse-Wishart prior
#' distribution on the variance-covariance matrix of the process error (W).
#' 
#' df_W: (numeric scalar) The degees of freedom for the inverse-Wishart prior
#' distribution on the variance-covariance matrix of the process error (W).
#' 
#' L: (numeric scalar) The period of the Fourier series approximation.
#' The spatial locations and knot locations are rescaled
#' to range from -L/4 to L/4 because the Fourier decomposition assumes that
#' the spatial surface is periodic.
#' 
#' smoothing: (numeric scalar) Controls the degree of smoothing in the 
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
#' @examples
#' # Create example data
#' num_time_points <- 5
#' spatial_locations <- expand.grid(seq(10), seq(10))
#' num_spatial_locations <- nrow(spatial_locations)
#' z <- rnorm(num_time_points * num_spatial_locations)
#' Y <- matrix(z, nrow=num_spatial_locations, ncol=num_time_points)
#' 
#' # Basic IDE model with one kernel
#' mod <- dstm_ide(Y, spatial_locations)
#' predict(mod)
#' summary(mod)
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
dstm_ide <- function(Y, locs=NULL, knot_locs=NULL, proc_error = "IW", J=3L,
                     n_samples = 10L, sample_sigma2 = TRUE,
                     verbose = FALSE, params = NULL) {
  
  # Error checking for J, L, locs, knot_locs
  if (is.null(locs)) stop("locs must be specified for IDE models")
  if (class(locs) == "data.frame") locs <- as.matrix(locs)
  if (class(locs) != "matrix") stop("locs must be data.frame or matrix")
  if (ncol(locs) != 2) stop("locs must have 2 columns")

  J <- as.integer(J)
  if (is.null(J) || is.na(J) || J<1) stop("J must be an integer > 0")
  P <- (2*J + 1)^2

  L <- params[["L"]] %else% 4
  check.numeric.scalar(L)
  
  smoothing <- params[["smoothing"]] %else% 1
  check.numeric.scalar(smoothing)
  
  SV <- is.numeric(knot_locs)
  if (!(SV || is.null(knot_locs))) stop("knot_locs must be numeric or null")
  
  # Center spatial locations to have range of L/2
  locs <- center_all(locs, L/2)
  
  # Process Model; creates kernel parameters
  locs_dim <- ncol(locs)
  
  # kernel_samples_per_iter
  kernel_samples_per_iter <- params[["kernel_samples_per_iter"]] %else% 1
  check.numeric.scalar(kernel_samples_per_iter, x_min=1, strict_inequality=FALSE)
  kernel_samples_per_iter <- as.integer(kernel_samples_per_iter)
  
  # mu_kernel_mean
  mu_kernel_mean <- params[["mu_kernel_mean"]] %else% rep(0, locs_dim)
  check.numeric.vector(mu_kernel_mean, locs_dim, dim_name="ncol(locs)")
  mu_kernel_mean <- as.matrix(mu_kernel_mean)
  
  # mu_kernel_var
  mu_kernel_var <- params[["mu_kernel_var"]] %else% diag(L/4, locs_dim)
  check.cov.matrix(mu_kernel_var, locs_dim, dim_name="ncol(locs)")
  
  # proposal_factor_mu
  proposal_factor_mu <- params[["proposal_factor_mu"]] %else% ifelse(SV, 2/5, 1)
  check.numeric.scalar(proposal_factor_mu)
  
  # Sigma_kernel_df
  k <- 10
  Sigma_kernel_df <- params[["Sigma_kernel_df"]] %else% (k * locs_dim)
  check.numeric.scalar(Sigma_kernel_df, x_min=locs_dim-1)
  
  # Sigma_kernel_scale
  Sigma_kernel_scale <- params[["Sigma_kernel_scale"]] %else% 
                          diag(Sigma_kernel_df*L/20, locs_dim)
  check.cov.matrix(Sigma_kernel_scale, locs_dim, dim_name="ncol(locs)")
  
  # proposal_factor_Sigma
  proposal_factor_Sigma <- params[["proposal_factor_Sigma"]] %else%
                           ifelse(SV, 1/12, 1)
  check.numeric.scalar(proposal_factor_Sigma)
  
  # Error checking for knot_locs
  # And adjustment to above quantities in spatially varying case
  n_knots <- 1
  K  <- array(0, dim=c(1, 1, 1))
  
  if (SV) {
    # prepare knot_locs
    if (length(knot_locs) > 1) knot_locs <- center_all(knot_locs, L/2)
    else knot_locs <- gen_grid(as.integer(knot_locs), L/2)
    
    # Modify kernel parameters
    n_knots <- nrow(knot_locs)
    mu_kernel_mean <- rep(1, n_knots) %x% matrix(mu_kernel_mean, nrow=1)
    mu_kernel_var <- mu_kernel_var %x% diag(n_knots)
    
    # Create K matrix
    K <- pdist::pdist(knot_locs, locs)
    K <- as.matrix(K)
    K <- exp(-K / smoothing)
    K <- apply(K, 2, function(x) x / sum(x))
    K <- t(K) # Makes K of dimension n_locs by n_knots
    K <- array(K, dim=c(dim(K), 1))
    
  }
  
  Sigma_kernel_scale <- array(1, dim=c(1, 1, n_knots)) %x%
                        Sigma_kernel_scale
  
  
  # Process remaining params
  new_params <- process_common_params(params, proc_error, P, sample_sigma2)
  scalar_params <- c(alpha_lambda=new_params[["alpha_lambda"]], 
                     beta_lambda=new_params[["beta_lambda"]],
                     alpha_sigma2=new_params[["alpha_sigma2"]], 
                     beta_sigma2=new_params[["beta_sigma2"]],
                     sigma2=new_params[["sigma2"]], 
                     df_W=new_params[["df_W"]], 
                     kernel_samples_per_iter=kernel_samples_per_iter,
                     proposal_factor_mu=proposal_factor_mu,
                     proposal_factor_Sigma=proposal_factor_Sigma,
                     Sigma_kernel_df=Sigma_kernel_df,
                     SV=SV, J=J, L=L)

  # Fit the model
  results <- ide(Y, locs, new_params[["m_0"]], new_params[["C_0"]],
                 mu_kernel_mean, mu_kernel_var, K, Sigma_kernel_scale, 
                 new_params[["C_W"]], scalar_params, n_samples, verbose)
  
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
  kernel_updates <- n_samples * kernel_samples_per_iter
  results[["mu_acceptance_rate"]] <- results[["mu_acceptances"]] / kernel_updates
  results[["Sigma_acceptance_rate"]] <- results[["Sigma_acceptances"]] / kernel_updates
  results[["mu_acceptances"]] <- NULL
  results[["Sigma_acceptances"]] <- NULL
  if (SV) results[["knot_locs"]] <- knot_locs
  
  class(results) <- c("dstm_ide" , "dstm", "list")
  attr(results, "proc_model") <- "ide"
  attr(results, "proc_error") <- proc_error
  attr(results, "sample_sigma2") <- sample_sigma2

  return(results)
}

.onUnload <- function (libpath) {
  library.dynam.unload("ideq", libpath)
}