# Methods for dstm objects

predict.dstm <- function(x, K = 1, only_K = FALSE, return_ys = TRUE,
                         return_thetas = FALSE, burnin = NULL) {
  if (is.null(burnin)) {
    if (is.null(x$burnin)) {
      message("INFO: Using all samples except starting values")
      burnin <- 1
    } else {
      burnin <- x$burnin
      message(paste("INFO: Not using first x$burnin =", burnin,
                    "samples (including starting values)"))
    }
  } else {
    message(paste("INFO: Not using first" , burnin,
                  "samples (including starting values)"))
  }

  # NOTE: burnin is 1 by default to remove starting values
  if (burnin < 1) stop("burnin must be greater than 0")
  S  <- nrow(x[["F"]])
  p  <- dim(x[["theta"]])[1]
  Tp1 <- dim(x[["theta"]])[2] # T plus 1
  n_samples <- dim(x[["theta"]])[3] - burnin
  idx <- seq(n_samples) + burnin

  # Calculate thetas
  thetas <- array(NA, dim = c(p, n_samples, K))

  # Functions to generate random process error
  get_W_chol <- function(W) lapply(W, function(M) t(chol(M)))
  get_w <- function(W_chol_) sapply(W_chol_, function(M) M %*% rnorm(p))

  # Calculate W
  if (attr(x, "proc_error") == "IW") {
    W <- lapply(seq(n_samples), function(i) chol2inv(x[["W_inv"]][,,i+burnin]))
  }
  else if (attr(x, "proc_error") == "discount") {
    if (attr(x, "proc_model") %in% c("AR", "Full")) {
      get_W <- function(i) {
        x[["lambda"]][i] *
          x[["G"]][,,i] %*% x[["C_T"]][,,i] %*% t(x[["G"]][,,i])
      W <- lapply(idx, get_W)
      }
    }
    else { # RW process model
      get_W <- function(i) x[["lambda"]][i] * x[["C_T"]][,,i]
      W <- lapply(idx, get_W)
    }
  }
  else {
    stop(paste("predict.dstm is not implemented for attr(x, \"proc_model\") == \"",
               attr(x, "proc_error"), "\""))
  }

  # Generate random process error
  W_chol <- get_W_chol(W)
  w <- get_w(W_chol)

  # Sample thetas for time T+1
  if (attr(x, "proc_model") %in% c("AR", "Full")) {
    get_thetas <- function(i) {
      j <- i+burnin
      x[["G"]][,,j] %*%
        x[["theta"]][,Tp1,j] + w[,i]
    }
    thetas[,,1] <- sapply(seq(n_samples), get_thetas)

  } else { # RW process model
    thetas[,,1] <- x[["theta"]][,Tp1,idx] + w
  }

  # Sample for T+2 up to T+K
  if (K > 1) {
    for (k in seq(2, K)) {
      # Need to recalculate W if using discount factors
      if (attr(x, "proc_error") == "discount") {
        W <- lapply(seq_along(W), function(i) x[["lambda"]][i+burnin] * W[[i]])
        W_chol <- get_W_chol(W)
      }

      # New process error
      w <- get_w(W_chol)

      # Generate thetas for T+k
      if (attr(x, "proc_model") %in% c("AR", "Full")) {
        get_thetas <- function(i) {
          j <- i + burnin
          x[["G"]][,, j] %*% thetas[,i,k-1] + w[, i]
        }
        thetas[,,k] <- sapply(seq(n_samples), get_thetas)
      } else { # RW process model
        thetas[,,k] <- thetas[,,k-1] + w
      }
    }
  }

  # Calculate ys
  if (return_ys) {
    # Calculate standard deviation for observation model
    if (attr(x, "sample_sigma2")) {
      my_sd <- rep(sqrt(x[["sigma2"]][idx]), each = p)
    }
    else {
      my_sd <- sqrt(x[["sigma2"]])
    }

    # Function to get predicted y values for a given time period
    get_preds <- function(k) x[["F"]] %*% thetas[,,k] + rnorm(n_samples * S, 0, my_sd)

    # Get predicted y values for requested time period
    if (only_K || K < 2) {
      ys <- get_preds(K)
    }
    # Get predicted y values for all time period <= T+K
    else {
      ys <- array(NA, dim = c(S,n_samples,K))
      for (k in seq(K)) {
        ys[,,k] <- get_preds(k)
      }
    }

    # Create output list depending on whether user wants thetas
    if (return_thetas) {
      if (only_K || K < 2) thetas <- thetas[,,K]
      results <- list(ys = ys, thetas = thetas)
    } else {
      results <- ys
    }

  # Create output for case when user does not want ys
  } else {
    if (only_K || K < 2) {
    results <- thetas[,,K]
    }
    else {
    results <- thetas
    }
  }

  return(results)
}

summary.dstm <- function(x, object_name = deparse(substitute(x))) {
  cat("Summary for dstm object \`", object_name, "\`\n", sep = "")
  cat("Observation model: `", attr(x, "obs_model"), "`\n", sep="")
  cat("Process model: `", attr(x, "proc_model"), "`\n", sep="")
  cat("Process error: `", attr(x, "proc_error"), "`\n", sep="")
  if (attr(x, "sample_sigma2")) {
    cat("sigma2 was sampled\n\n")
  } else {
    cat("sigma2 was taken as a constant\n\n")
  }

  cat("List elements (in order) are as follows:\n")
  cat(names(x), "\n")

  # Vector Objects
  cat("\nVector Objects:\n")
  numeric_idx <- which(sapply(x, function(y) is.vector(y)))
  numeric_summary <- matrix(NA, nrow = length(numeric_idx), ncol = 8)
  my_probs = c(0.0, 0.25, 0.5, 0.75, 1.0)
  counter <- 1
  for (i in numeric_idx) {
    numeric_summary[counter, 1]   <- length(x[[i]])
    numeric_summary[counter, 2]   <- mean(x[[i]])
    numeric_summary[counter, 3]   <- var(x[[i]])
    numeric_summary[counter, 4:8] <- quantile(x[[i]], probs = my_probs)
    counter <- counter + 1
  }
  colnames(numeric_summary) <- c("Length", "Mean", "Var",
                               paste0(as.character(my_probs * 100), "%"))
  rownames(numeric_summary) <- names(x)[numeric_idx]
  print(numeric_summary)

  # Other objects
  cat("\n\nOther Objects:\n")
  other_idx <- setdiff(seq_along(x), numeric_idx)
  other_summary <- matrix(NA, nrow = length(other_idx), ncol = 4)
  counter <- 1
  for (i in other_idx) {
    other_summary[counter, 1] <- class(x[[i]])
    if (is.matrix(x[[i]]) || is.array(x[[i]])) {
      dims_i <- dim(x[[i]])
      other_summary[counter, 2:(1 + length(dims_i))] <- dims_i
    }
    counter <- counter + 1
  }
  colnames(other_summary) <- c("class", "dim 1", "dim 2", "dim 3")
  rownames(other_summary) <- names(x)[other_idx]
  print(other_summary)
}

# print.dstm is just a wrapper for summary.dstm
print.dstm <- function(x, display = deparse(substitute(x))) {
  summary.dstm(x, object_name = display)
}
