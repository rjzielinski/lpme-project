weight_seq <- function(x.obs, mu, sigma, epsilon = 0.001, max.iter = 1000) {
  # x.obs is the dataset of interest
  # there are n observations, and each observation is a D-dimensional point
  # x.obs is a n-by-D matrix

  # mu is a vector of the knots in a mixture density estimation
  # sigma is the bandwidth of the density estimation
  # epsilon is a predetermined tolerance of the Euclidean distance between thetas in two consecutive steps
  # max.iter is a predetermined upper bound on the number of steps

  n.D <- dim(x.obs)
  n <- n.D[1]
  D <- n.D[2]
  N <- dim(mu)[1]

  A <- matrix(NA, nrow = n, ncol = N)
  for (j in 1:N) {
    A.prepare <- function(x) {
      return(ker(x, mu[j, ], sigma))
    }
    A[, j] <- apply(x.obs, 1, A.prepare)
  }

  theta.old <- rep(1 / N, N)
  abs.diff <- 10 * epsilon
  count <- 0
  lambda.hat.old <- c(n, rep(-1, D))

  while((abs.diff > epsilon) & (count <= max.iter)) {
    W <- t(t(A) * theta.old)
    W <- W / apply(W, 1, sum)
    w <- apply(W, 2, sum)

    flambda <- function(lambda) {

      denom.temp <- apply(
        t(t(cbind(rep(1, dim(mu)[1]), mu)) * lambda),
        1,
        sum
      )
      num.temp <- mu * w

      f1 <- sum(w / denom.temp)
      f2 <- apply(
        num.temp * (as.vector(1 / denom.temp)),
        2,
        sum
      )
      f <- dist_euclidean(f1, 1) + dist_euclidean(f2, apply(x.obs, 2, mean))
      return(f)
    }

    lambda.hat <- nlm(flambda, lambda.hat.old, iterlim = 1000)$estimate

    theta.new <- w /
      apply(
        t(t(cbind(rep(1, dim(mu)[1]), mu)) * lambda.hat),
        1,
        sum
      )

    abs.diff <- dist_euclidean(theta.new, theta.old)
    if (is.na(abs.diff)) {
      abs.diff <- 0
      theta.new <- theta.old
    }

    theta.old <- pmax(theta.new, 0)
    theta.old <- pmin(theta.old, 1)
    count <- count + 1
    lambda.hat.old <- lambda.hat
  }

  theta.hat <- pmax(theta.new, 0)
  theta.hat <- pmin(theta.hat, 1)
  return(theta.hat)
}