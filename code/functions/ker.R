ker <- function(x, mu, sigma) {
  # Smoothing kernel for density estimation
  yseq <- dnorm((x - mu) / sigma)
  return((sigma ^ {-length(x)}) * prod(yseq))
}