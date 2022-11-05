hdmde <- function(x.obs, N0, alpha, max.comp) {

  # "x.obs" is the data set of interest.
  # There are n observations, and each observation is a D-dimensional point.
  # x.obs is a n-by-D matrix.
  # "N0" is a predetermined lower bound for N - the number of density components
  # "alpha" is the predetermined confidence level.
  # "max.comp" is the upper bound of the number of components in the desired mixture density.

  zalpha <- qnorm(1 - alpha / 2)
  n.D <- dim(x.obs) # "x.obs" is a n by D matrix. Each row of it denotes a data point in R^D.
  n <- n.D[1]
  D <- n.D[2]
  N <- N0

  km <- kmeans(x.obs, N, iter.max = 100, nstart = 100) # Use k-means clustering to get N clusters of x.obs.
  mu <- km$centers # Set the centers of the N clusters as the means of the density components.

  sigma.vec <- rep(NA, N) # The following block estimates \sigma_N.
  for(j in 1:N) {
    index.temp <- which(km$cluster == j)
    xi.j <- matrix(x.obs[index.temp, ], nrow = length(index.temp))
    sig.prepare <- function(x) {
      return((dist_euclidean(x, mu[j, ])) ^ 2)
    }
    s <- apply(xi.j, 1, sig.prepare)
    sigma.vec[j] <- mean(s)
  }
  sig <- sqrt(mean(sigma.vec) / dim(x.obs)[2])

  # It gives an estimation of theta_j's with weight.seq().
  theta.hat <- weight_seq(x.obs, mu, sig)

  # The following block gives an approximation to the underlying density function of interest.
  # This estimation is of the form of weights times scaled kernels.
  f.test <- function(x) {
    fun.prepare <- function(t) {
      return(ker(x, t, sig))
    }
    comp.vec <- apply(mu, 1, fun.prepare)
    return(sum(theta.hat * comp.vec))
  }

  p.old <- apply(x.obs,1,f.test) # The second/negative term p_N in Delta.hat.

  test.rejection <- 1
  while ((test.rejection == 1) & (N <= min(n, max.comp))) {

    if (max.comp > N0) {
      N <- N + 1
    }

    ##################################################
    # The following is a repetition of the codes above.
    km <- kmeans(x.obs, N, iter.max = 100, nstart = 100)
    mu <- km$centers

    sigma.vec <- rep(NA, N)
    for(j in 1:N) {
      index.temp <- which(km$cluster == j)
      xi.j <- matrix(x.obs[index.temp, ], nrow = length(index.temp))
      sig.prepare <- function(x) {
        return(dist_euclidean(x, mu[j, ]) ^ 2)
      }
      s <- apply(xi.j, 1, sig.prepare)
      sigma.vec[j] <- mean(s)
    }
    sig <- sqrt(mean(sigma.vec) / dim(x.obs)[2])

    theta.hat <- weight_seq(x.obs, mu, sig)

    f.test <- function(x) {
      fun.prepare <- function(t) {
        return(ker(x, t, sig))
      }
      comp.vec <- apply(mu, 1, fun.prepare)
      return(sum(theta.hat * comp.vec))
    }

    # The part above is a repetition.
    ##################################################

    p.new <- apply(x.obs, 1, f.test) # The first/positive term p_{N+1} in Delta.hat.
    delta.hat <- p.new - p.old # Delta.hat
    sigma.hat.sq <- mean((delta.hat - mean(delta.hat)) ^ 2)
    Z.I.N <- sqrt(dim(x.obs)[1]) * mean(delta.hat) / sqrt(sigma.hat.sq)

    if (
      (Z.I.N <= zalpha) &
      (Z.I.N >= -zalpha) &
      (!is.na(Z.I.N))
    ) {
      test.rejection=0
    }
    p.old <- p.new
  }

  f <- f.test

  resp <- list(
    estimating.pdf = f,
    theta.hat = theta.hat,
    mu = mu,
    N = N,
    k.means.result = km,
    sigma = sig,
    Z.I.N = Z.I.N
  )

  # Output:
  # "estimating.pdf" is the derived mixture density approximating an underlying desnity.
  # "theta.hat" is a vector of weights for knots of this mixture density.
  # "mu" is a vector of knots of this mixture density.
  # "N" is the number of knots of this mixture density.
  # "sigma" is the variance shared by the components of this mixture density.
  # "Z.I.N" is the statistic determining the size of N.
  # "k.means.result" gives the result list of "kmeans(obs.x,N)".

  return(resp)
}
