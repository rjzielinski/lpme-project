ker <- function(x, mu, sigma) {
  # Smoothing kernel for density estimation
  #yseq <- sapply((x - mu) / sigma, dnorm)
  yseq <- dnorm((x - mu) / sigma)
  return((sigma ^ {-length(x)}) * prod(yseq))
}