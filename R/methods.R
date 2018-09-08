# Methods for dstm objects

predict.dstm <- function(x, K = 1, only_K = FALSE, return_ys = TRUE,
                         return_thetas = FALSE) {
  S  <- nrow(x[["F"]])
  p  <- dim(x[["theta"]])[1]
  Tp1 <- dim(x[["theta"]])[2] # T plus 1
  n_samples <- dim(x[["theta"]])[3]

  # Calculate thetas
  thetas <- array(NA, dim = c(p, K, n_samples))

  # proc_error = IW
  if (attr(x, "proc_error") == "IW") {
    W_chol <- lapply(1:n_samples + 1, function(i) solve(chol(x[["W_inv"]][ , , i])))
    w <- sapply(1:n_samples, function(i) W_chol[[i]] %*% rnorm(p))

    if (attr(x, "proc_model") %in% c("AR", "Full")) {
      thetas[, 1,] <- sapply(1:n_samples, function(i) x[["G"]][,, i + 1] %*%
                                            x[["theta"]][, Tp1, i] + w[, i])
    } else {
        thetas[, 1,] <- x[["theta"]][ , Tp1, ] + w
    }

    if (K > 1) {
      for (k in 2:K) {
        if (attr(x, "proc_model") %in% c("AR", "Full")) {
          thetas[, k,] <- sapply(1:n_samples, function(i) x[["G"]][,, i + 1] %*%
                                              thetas[, k - 1, i] + w[, i])
        } else {
            thetas[, k,] <- thetas[, k - 1,] + w
        }
      }
    }

  # proc_error != IW
  } else {
    stop("predict.dstm only implemented for attr(x, \"proc_model\") == \"IW\"")
  }

  # Calculate ys
  if (return_ys) {
    if (attr(x, "sample_sigma2")) {
      my_sd <- rep(sqrt(x[["sigma2"]][-1]), each = p)
    } else {
      my_sd <- sqrt(x[["sigma2"]])
    }

    if (only_K) {
      ys <- x[["F"]] %*% thetas[, K,] + rnorm(n_samples * S, 0, my_sd)
    } else {
      ys <- sapply(1:K, function(k) x[["F"]] %*% thetas[, k,] +
                                    rnorm(n_samples * S, 0, my_sd))
      if (K > 1) {
        ys <- array(ys, dim = c(S, n_samples, K))
      } else {
        ys <- matrix(ys, nrow = S, ncol = n_samples)
      }
    }

    if (return_thetas) {
      if (K < 2) thetas <- thetas[, 1,]
      results <- list(ys = ys, thetas = thetas)
    } else {
      results <- ys
    }

  # Only thetas
  } else if (only_K) {
    results <- thetas[, K, ]
  } else {
    if (K < 2) thetas <- thetas[, 1,]
    results <- thetas
  }

  return(results)
}

summary.dstm <- function(x, object_name = deparse(substitute(x))) {
  cat("Summary for dstm object \`", object_name, "\`\n\n", sep = "")
  cat("List elements (in order) are as follows:\n")
  print(names(x))

  # Vector Objects
  cat("\nVector Objects:\n")
  numeric_idx <- which(sapply(x, function(y) length(dim(y)) < 3 && ncol(y) == 1))
  numeric_summary <- matrix(NA, nrow = length(numeric_idx), ncol = 7)
  my_probs = c(0.0, 0.025, 0.5, 0.975, 1.0)
  counter <- 1
  for (i in numeric_idx) {
    numeric_summary[counter, 1]   <- mean(x[[i]])
    numeric_summary[counter, 2]   <- var(x[[i]])
    numeric_summary[counter, 3:7] <- quantile(x[[i]], probs = my_probs)
    counter <- counter + 1
  }
  colnames(numeric_summary) <- c("Mean", "Var",
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
