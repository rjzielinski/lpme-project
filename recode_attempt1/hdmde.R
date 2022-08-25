hdmde2 <- function(x.obs, N0, alpha, max.comp, p.old = NULL) {
  
  zalpha <- qnorm(1 - alpha / 2)
  n.D <- dim(x.obs)
  n <- n.D[1]
  D <- n.D[2]
  N <- N0
  
  print(N)
  
  km <- kmeans(
    x.obs,
    N,
    nstart = 100
  )  
  mu <- km$centers
  
  # calc_sigma <- function(km, j) {
  #   index.temp <- which(km$cluster == j)
  #   xi.j <- x.obs[index.temp, ]
  #   sig.prepare <- function(x) {
  #     return((dist.euclidean(x, mu[j, ]))^2)
  #   }
  #   s <- apply(xi.j, 1, sig.prepare)
  #   return(mean(s))
  # }
  
  # sigma.vec <- sapply(1:N, calc_sigma, km = km)
  
  sigma.vec <- rep(NA, N)
  # for (idx in 1:length(sort(unique(km$cluster)))) {
  #   sigma.vec[idx] <- calc_sigma(km, sort(unique(km$cluster))[idx])
  # }
  for (j in 1:N) {
    index.temp <- which(km$cluster == j)
    xi.j <- x.obs[index.temp, ]
    sig.prepare <- function(x) {
      return((dist.euclidean(x, mu[j, ]))^2)
    }
    s <- apply(xi.j, 1, sig.prepare)
    sigma.vec[j] <- mean(s)
  }
  
  sig <- sqrt(mean(sigma.vec) / dim(x.obs)[2])
  
  theta.hat <- weight.seq(x.obs, mu, sig)
  
  f.test <- function(x) {
    comp.vec <- ker(x, mu, sig)
    return(sum(theta.hat * comp.vec))
  }
  
  p.old <- apply(x.obs, 1, f.test)
  
  if (is.null(p.old)) {
    p.new <- apply(x.obs, 1, f.test)
    test.rejection <- 1
  } else {
    p.new <- apply(x.obs, 1, f.test)
    delta.hat <- p.new - p.old
    sigma.hat.sq <- mean((delta.hat - mean(delta.hat))^2)
    Z.I.N <- sqrt(dim(x.obs)[1]) * mean(delta.hat) / sqrt(sigma.hat.sq)
    if ((Z.I.N <= zalpha) & (Z.I.N >= -zalpha) & (!is.na(Z.I.N))) {
      test.rejection <- 0
    } else {
      test.rejection <- 1
    }
  }
  
  if (test.rejection == 0 | N > min(n, max.comp)) {
    resp <- list(
      estimating.pdf = f.test,
      theta.hat = theta.hat,
      mu = mu,
      N = N,
      k.means.result = km,
      sigma = sig,
      Z.I.N = Z.I.N
    )
    
    # Output:
    # estimating.pdf is the derived mixture density approximating an underlying density
    # theta.hat is a vector of weights for knots of this mixture density
    # mu is a vector of knots of this mixture density
    # N is the number of knots of this mixture density
    # sigma is the variance shared by the components of this mixture density
    # Z.I.N is the statistic determining the size of N
    # k.means.result gives the result list of kmeans(obs.x, N)
    
    return(resp)
  } else {
    return(hdmde2(x.obs, N + 1, alpha, max.comp, p.new))
  }
}
