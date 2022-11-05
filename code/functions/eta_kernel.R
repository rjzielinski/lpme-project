eta_kernel <- function(t, lambda) {
  # reproducing kernels associated with Sobolev space D^{-2}L^2(R^d)
  norm_val <- norm_euclidean(t)
  if (lambda %% 2 == 0) {
    if (norm_val == 0) {
      y <- 0
    } else {
      y <- (norm_val ^ lambda) * log(norm_val)
    }
  } else {
    y <- norm_val ^ lambda
  }
  return(y)
}