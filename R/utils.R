# Functions for dstm_ide()
center_col <- function(x, L) {
  x_range <- diff(range(x))
  x <- L/x_range * x
  x <- x + (L/2 - max(x))
  x
}

center_all <- function(x, L) {
  for ( i in seq(ncol(x)) ) {
    x[, i] <- center_col(x[, i], L)
  }
  x
}

gen_grid <- function(kernel_locs, L) {
  increment <- L / (kernel_locs)
  extreme <- (L - increment) / 2
  x <- seq(-extreme, extreme, increment)
  g <- expand.grid(x, x)
  g <- as.matrix(g)
  colnames(g) <- c("x", "y")
  g
}

# Functions for predict.dstm()
update_C <- function(C_prev, G) {
  C_new <- array(-1, dim=dim(C_prev))
  n_samples <- dim(C_new)[3]
  for (i in seq(n_samples)) C_new[,,i] <- G[,,i] %*% C_prev[,,i] %*% t(G[,,i])
  C_new
}

calc_W <- function(lambda, C_T) {
  W <- array(-1, dim=dim(C_T))
  n_samples <- dim(C_T)[3]
  for (i in seq(n_samples)) W[,,i] <- lambda[i] * C_T[,,i]
  W
}

next_thetas <- function(thetas_prev, G, W) {
  n_samples <- ncol(thetas_prev)
  E_t <- sapply(seq(n_samples), function(i) G[,,i] %*% thetas_prev[,i])
  sapply(seq(n_samples), function(i) mvtnorm::rmvnorm(1, E_t[,i], W[,,i]))
}
