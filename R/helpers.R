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
