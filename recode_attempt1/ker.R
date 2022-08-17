ker <- function(x, mu, sigma) {
  # Smoothing kernel for density estimation
  # We applied Gaussian kernel
  require(Rfast)
  yseq <- dnorm((x - t(mu)) / sigma)
  return((sigma ^ {-length(x)}) * colprods(yseq))
}