# Methods for dstm objects

#' Predict Method for DSTM Fits
#' @description 
#' Generates samples from the posterior predictive distribution 
#' at future time points for 
#' (1) the observation vector and (2) the state vector
#' @param x
#' A `dstm` object
#' 
#' @param K 
#' (integer scalar) The number of future time periods for which to generate
#' predictions
#' 
#' @param only_K 
#' (logical scalar) Whether to return predictions for time period T+K only
#' (as opposed to T+1, T+2, ..., T+K)
#' 
#' @param return_ys
#' (logical scalar) Whether to return samples from the posterior predictive
#' distribution of the observation vector (ys)
#' 
#' @param return_thetas
#' (logical scalar) Whether to return samples from the posterior predictive
#' distribution of the state vector (thetas)
#' 
#' @param burnin
#' (integer scalar) The number of samples to discard as burn-in.
#' If x$burnin exists, this argument will override it.
#' 
#' @details
#' The posterior predictive samples are returned in a matrix or 3-D array,
#' depending on whether samples from multiple time points are requested.
#' The dimensions are always in the following order:
#' 
#' 1. p: The index of each parameter within the state vector.
#' 
#' 2. t: The time index
#' 
#' 3. n_samples: The sample number
#' @export
#' @examples
#' num_time_points <- 5
#' spatial_locations <- expand.grid(seq(10), seq(10))
#' num_spatial_locations <- nrow(spatial_locations)
#' z <- rnorm(num_time_points * num_spatial_locations)
#' Y <- matrix(z, nrow=num_spatial_locations, ncol=num_time_points)
#' 
#' # ide example
#' mod_ide <- dstm_ide(Y, spatial_locations)
#' predict(mod_ide)
#' predict(mod_ide, K=4, return_thetas=TRUE)
#' 
#' # eof example
#' mod_eof <- dstm_eof(Y)
#' predict(mod_eof, K=2, only_K=TRUE, burnin=5)
predict.dstm <- function(x, K = 1, only_K = FALSE, return_ys = TRUE,
                         return_thetas = FALSE, burnin = NULL) {
  # Argument burnin is used if provided
  # If not, we use x$burnin if available
  # Otherwise, we use no burnin
  if (is.null(burnin)) {
    burnin <- if ( is.null(x$burnin) ) 0 else x$burnin
  }
  if (burnin < 0) stop("burnin must be non-negative")
  if (!return_ys && !return_thetas)
    stop("return_ys and/or return_thetas must be true")
  
  # Model meta information
  Discount <- attr(x, "proc_error") == "Discount"
  RW <- attr(x, "proc_model") == "RW"
  
  # Dimensions
  S  <- nrow(x[["F"]])
  P  <- dim(x[["theta"]])[1]
  Tp1 <- dim(x[["theta"]])[2] # T plus 1
  n_samples <- dim(x[["theta"]])[3] - burnin
  idx_post_burnin <- seq(n_samples) + burnin
  
  # Create copies of objects needed for sampling
  thetas_prev <- x[["theta"]][,Tp1,idx_post_burnin]
  sigma2 <- x[["sigma2"]][idx_post_burnin]
  if (RW) {
    G <- array(1, dim=c(1, 1, n_samples)) %x% diag(P)
  } else {
    G <- x[["G"]][,,idx_post_burnin]
  }
  
  if (Discount) {
    lambda <- x[["lambda"]][idx_post_burnin]
    C_T <- x[["C_T"]][,,idx_post_burnin]
  } else {
    W <- x[["W"]][,,idx_post_burnin]
  }
  
  # Step 1: Sample thetas from posterior predictive distribution
  thetas <- array(NA, dim = c(P, K, n_samples))
  
  # Create W for discount models
  if (Discount) {
    C_Tpk <- update_C(C_T, G)
    W <- calc_W(lambda, C_Tpk)
  }
  
  thetas[,1,] <- next_thetas(thetas_prev, G, W)
  if (K > 1) {
    for (k in seq(2, K)) {
      if (Discount) {
        C_Tpk <- update_C(C_Tpk, G)
        W <- calc_W(lambda, C_Tpk)
      }
      thetas[,k,] <- next_thetas(thetas[,k-1,], G, W)
    }
  }

  # Calculate ys
  if (return_ys) {
    # Calculate standard deviation for observation model
    sample_sigma2 <- is.logical(attr(x, "sample_sigma2")) &&
                     attr(x, "sample_sigma2")
    if (sample_sigma2) my_sd <- rep(sqrt(sigma2), each=P)
    else my_sd <- sqrt(x[["sigma2"]])

    # Function to get predicted y values for a given time period
    get_preds <- function(k) x[["F"]] %*% thetas[,k,] + rnorm(n_samples*S, sd=my_sd)

    # Get predicted y values for requested time period
    if (only_K || K < 2) {
      ys <- get_preds(K)
    }
    else { # Get predicted y values for all time period <= T+K
      ys <- array(NA, dim = c(S, n_samples, K))
      for (k in seq(K)) ys[,,k] <- get_preds(k)
    }

    # Create output list depending on whether user wants thetas
    if (return_thetas) {
      if (only_K || K < 2) thetas <- thetas[,K,]
      results <- list(ys = ys, thetas = thetas)
    } else results <- ys

  # Create output for case when user does not want ys
  } else {
    if (only_K || K < 2) results <- thetas[,K,]
    else results <- thetas
  }

  return(results)
}

#' Summary Method for DSTM Fits
#' @description 
#' Prints summary information for `dstm` objects
#' @param x A `dstm` object
#' @param object_name The name to be printed in the summary (if desired)
#' @examples
#' num_time_points <- 5
#' num_spatial_locations <- 100
#' z <- rnorm(num_time_points * num_spatial_locations)
#' Y <- matrix(z, nrow=num_spatial_locations, ncol=num_time_points)
#' rw_model <- dstm_eof(Y, proc_model="RW") # Random walk model
#' summary(rw_model)
#' 
#' # Identically
#' print(rw_model)
#' rw_model
#' @export
summary.dstm <- function(x, object_name = deparse(substitute(x))) {
  cat("Summary for dstm object \`", object_name, "\`\n", sep = "")
  cat("Process model: `", attr(x, "proc_model"), "`\n", sep="")
  cat("Process error: `", attr(x, "proc_error"), "`\n", sep="")
  if (attr(x, "sample_sigma2")) {
    cat("sigma2 was sampled\n\n")
  } else {
    cat("sigma2 was fixed\n\n")
  }

  cat("List elements (in order) are as follows:\n")
  cat(names(x), "\n")
  
  # Numeric Scalars
  numeric_bool <- sapply(x, function(y) is.vector(y) && is.numeric(y))
  scalar_bool <- sapply(x, function(y) length(y) == 1)
  scalar_idx <- which(numeric_bool &  scalar_bool)
  vector_idx <- which(numeric_bool & !scalar_bool)
  
  if ( length(scalar_idx) > 0 ) {
    cat("\nScalar Objects:\n")
    scalars <- numeric()
    for (i in scalar_idx) scalars <- c(scalars, x[[i]])
    names(scalars) <- names(x)[scalar_idx]
    print(scalars)
  }

  # Numeric Vectors
  if ( length(vector_idx) > 0 ) {
    cat("\nVector Objects:\n")
    vector_summary <- matrix(NA, nrow = length(vector_idx), ncol = 8)
    my_probs = c(0.0, 0.25, 0.5, 0.75, 1.0)
    counter <- 1
    for (i in vector_idx) {
      vector_summary[counter, 1]   <- length(x[[i]])
      vector_summary[counter, 2]   <- mean(x[[i]])
      vector_summary[counter, 3]   <- var(x[[i]])
      vector_summary[counter, 4:8] <- quantile(x[[i]], probs = my_probs)
      counter <- counter + 1
    }
    colnames(vector_summary) <- c("Length", "Mean", "Var",
                                 paste0(as.character(my_probs * 100), "%"))
    rownames(vector_summary) <- names(x)[vector_idx]
    print(vector_summary)
  }

  # Matrices/Arrays
  mat_arr_idx <- which(sapply(x, function(y) is.matrix(y) || is.array(y)))
  if ( length(mat_arr_idx) > 0 ) {
    cat("\n\nMatrices/Arrays:\n")
    mat_arr_summary <- matrix(NA, nrow = length(mat_arr_idx), ncol = 4)
    counter <- 1
    for (i in mat_arr_idx) {
      mat_arr_summary[counter, 1] <- class(x[[i]])
      dims_i <- dim(x[[i]])
      mat_arr_summary[counter, 2:(1 + length(dims_i))] <- dims_i
      counter <- counter + 1
    }
    colnames(mat_arr_summary) <- c("class", "dim 1", "dim 2", "dim 3")
    rownames(mat_arr_summary) <- names(x)[mat_arr_idx]
    mat_arr_summary <- as.data.frame(mat_arr_summary)
    print(mat_arr_summary)
  }
  
  # Lists
  list_idx <- which(sapply(x, function(y) is.list(y)))
  if ( length(list_idx) > 0 ) {
    cat("\n\nLists:\n")
    list_summary <- matrix(NA, nrow = length(list_idx), ncol = 4)
    counter <- 1
    for (i in list_idx) {
      list_summary[counter, 1] <- length(x[[i]])
      if (!grepl("param", names(x)[i])) {
        dims_i <- dim(x[[i]][[1]])
        list_summary[counter, 2:(1 + length(dims_i))] <- dims_i
      }
      counter <- counter + 1
    }
    colnames(list_summary) <- c("length", "dim 1", "dim 2", "dim 3")
    rownames(list_summary) <- names(x)[list_idx]
    print(list_summary)
  }
}

#' Print Method for DSTM Fits
#' @description 
#' Prints a summary for a `dstm` object by calling summary.dstm().
#' @param x A `dstm` object
#' @seealso summary.dstm
#' @export
print.dstm <- function(x, display = deparse(substitute(x))) {
  summary.dstm(x, object_name = display)
}
