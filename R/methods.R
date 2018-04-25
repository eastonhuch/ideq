# Methods for dstm objects

predict.dstm <- function(x) {
  print("not implemented. Ha!")
}
predict(dat_full)

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
summary(dat_full)

# print.dstm is just a wrapper for summary.dstm
print.dstm <- function(x, display = deparse(substitute(x))) {
  summary.dstm(x, object_name = display)
}
print(dat_full)
