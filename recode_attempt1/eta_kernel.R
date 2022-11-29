eta.kernel <- function(t, lambda) {
  norm_val <- norm.euclidean(t)
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