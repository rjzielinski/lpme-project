weight.seq <- function(x.obs, mu, sigma, epsilon=0.001, max.iter = 1000) {
  
  # x.obs is the dataset of interest.
  #       There are n observations, with each observation being a D-dimensional point.
  #       x.obs is an n x D matrix
  # mu is a vector of the knots in a mixture density estimation
  # sigma is the bandwidth of this density estimation
  # epsilon is a predetermined tolerance of the Euclidean distance between thetas in two consecutive steps
  # max.iter is a predetermined upper bound of the number of steps in this iteration
  
  require(Rfast)
  
  n.D <- dim(x.obs)
  n <- n.D[1]
  D <- n.D[2]
  N <- dim(mu)[1]
  
  A <- matrix(NA, ncol = N, nrow = n)
  
  for (i in 1:n) {
    A[i, ] <- ker(x.obs[i, ], mu = mu, sigma = sigma)
  }
  
  theta.old <- rep(1 / N, N) # The initial guess of the weights, theta_j
  abs.diff <- 10 * epsilon # Absolute value of the difference between "theta.new" and "theta.old"
  count <- 0 # Counting the number of steps in this iteration
  lambda.hat.old <- c(n, rep(-1, D)) # The initial guess of the Lagrangian multipliers
  
  while ((abs.diff > epsilon) & (count <= max.iter)) { # The iteration for computing desired thetas
    W <- t(t(A) * theta.old)
    W <- W / rowsums(W)
    w <- colsums(W)
    
    flambda <- function(lambda) {
      denom.temp <- rowsums(t(t(cbind(1, mu)) * lambda))
      num.temp <- mu * w
      
      f1 <- sum(w / denom.temp)
      f2 <- colsums(num.temp * as.vector(1 / denom.temp))
      f <- dist.euclidean(f1, 1) + dist.euclidean(f2, colmeans(x.obs))
      return(f)
    }
    
    lambda.hat <- nlm(flambda, lambda.hat.old, iterlim = 1000)$estimate
    
    # We set the Lagrangian multipliers in the previous step
    # as the initial guess in this step
    
    theta.new <- w / rowsums(t(t(cbind(1, mu)) * lambda.hat))
    abs.diff <- dist.euclidean(theta.new, theta.old)
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
