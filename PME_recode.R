#########################################################################
##########    Principal Manifold Estimation (PME) Algorithm   ###########
##########                  Kun (Michael) Meng                ###########
##########    Department of Biostatistics, Brown University   ###########
#########################################################################


########### Introduction ###########################################
####################################################################

# We propose an R-software function for d-dimensional estimating principle manifolds, 
# where d = 1, 2, 3, and manifolds embedded into D-dimensional Euclidean space with D 
# strictly larger than d.

# The proposed function is "PME()." To construct this function, we need some functions 
# for "high dimensional mixture density estimation (HDMDE)."

# This manuscript - whatever its name is - is organized as follows.
#   (i)   In Section 1, we define some basic functions.
#   (ii)  In Section 2, we construct a function for HDMDE. In this section, we first construct 
#         the "weight.seq" function. Then we apply this function to construct "hdmde" function.
#   (iii) Based on this "hdmde" function, "PME" estimation function is finally given in Section 3.
# The theoretical foundation of this R-software function is the following paper.
#
# K. Meng and A. Eloyan, Principal Manifolds: A framework Using Sobolev Spaces and Model
#                        Complexity Selection Using Mixture Densities

########### Section 0, Packages ####################################
####################################################################

#install.packages("MASS")
library(MASS)
#install.packages("Matrix")
library(Matrix)
library(memoise)
#install.packages("vegan")
library(vegan)
#install.packages("plot3D")
library(plot3D)
library(Rfast)
library(tidyverse)


########### Section 1, Some Basic Functions ########################
####################################################################

## Subsection 1.1, Functions for Euclidean metrics

# Norm function in an Euclidean space of any dimension
norm_euclidean <- function(x) {
  temp_vec <- as.numeric(x)
  norm_val <- sum(temp_vec ^ 2) %>% 
    sqrt()
  return(norm_val)
}

# Distance function in an Euclidean space of any dimension
dist_euclidean <- function(x, y) {
  return(norm_euclidean(x - y))
}

## Subsection 1.2, Kernels for minimization in a semi-normed space of Sobolev type

# Smoothing kernel for density estimation 
# (We applied Gaussian kernel.)
ker <- function(x, mu, sigma) {
  yseq <- sapply((x - mu) / sigma, dnorm)
  return((sigma ^ {-length(x)}) * prod(yseq))
}

# Reproducing Kernels associated with Sobolev space D^{-2}L^2(R^d)
eta_kernel <- function(t, lambda) {
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

## Subsection 1.3, Projection Index function
projection <- function(x, f, initial_guess) {
  DD <- function(t) {
    return(dist_euclidean(x, f(t)))
  }
  est <- nlm(DD, p = initial_guess)
  return(est$estimate)
}

# mem_projection <- memoise(projection)
mem_projection <- memoise(projection)

##### Section 2, High Dimensional Mixture Density Estimation #######
####################################################################

## Subsection 2.1 
# When \mu's and \sigma are given, the following function estimates \hat{\theta}'s.

weight_seq <- function(x_obs, mu, sigma, epsilon = 0.001, max_iter = 1000) { 
  
  # "x.obs" is the data set of interest. 
  #         There are n observations, and each observation is a D-dimensional point.
  #         x.obs is a n-by-D matrix.
  # "mu" is a vector of the knots in a mixture density estimation.
  # "sigma" is the bandwidth of this density estimation.
  # "epsilon" is a predetermined tolerance of the Euclidean distance between thetas in two consecutive steps.
  # "max.iter" is a predetermined upper bound of the number of steps in this iteration.
  
  n_D <- dim(x_obs)
  n <- n_D[1]
  D <- n_D[2]
  N <- dim(mu)[1]                     
  
  A <- matrix(NA, ncol = N, nrow = n)
  for(j in 1:N) {
    A_prepare <- function(x) {
      return(ker(x, mu[j, ], sigma))  
    }
    A[, j] <- apply(x_obs, 1, A_prepare) # A[i,j] is \psi_sigma (x_i-mu_j).
  }
  
  # The initial guess for weights, say \theta_j's.
  theta_old <- rep(1 / N, N)
  # Absolute value of the difference between "theta.new" and "theta.old".
  abs_diff <- 10 * epsilon 
  # Counting the number of steps in this iteration.
  count <- 0 
  # The initial guess of the Lagrangian multipliers
  lambda_hat_old <- c(n, rep(-1, D)) 
  
  # The iteration for computing desired theta's
  while((abs_diff > epsilon) & (count <= max_iter)) { 
    
    W <- t(t(A) * theta_old) # \theta_j^{(k)} \times \psi_\sigma(x_i-mu_j)
    W <- W / rowsums(W) # W[i,j] is the posterior probability of Z=j|X=x_i, say w_{i,j}(\theta.old).
    w <- colsums(W) # w[j] = \sum_{i=1}^n w_{i,j}
    
    flambda <- function(lambda) {
      # This function is for computing Lagrangian multipliers.
      denom_temp <- rowsums(
        t(t(cbind(rep(1, dim(mu)[1]), mu)) * lambda)
      )
      
      # \lambda_1+\lambda_2^T \mu_j, j=1,2,...,N.
      num_temp <- mu * w
      
      # \sum_{j=1}^N \frac{ w_ij }{ \lambda_1+\lambda_2^T \mu_j }
      f1 <- sum(w / denom_temp)
      f2 <- colsums(num_temp * as.vector(1 / denom_temp))
      f <- dist_euclidean(f1, 1) + dist_euclidean(f2, colmeans(x_obs))
      return(f) 
    }
    
    lambda_hat <- nlm(flambda, lambda_hat_old, iterlim=1000)$estimate              # The lagrangian multipliers.
    # We set the Lagrangian multipliers in the previous step 
    # as the initial guess in this step.
    theta_new <- w / rowsums(t(t(cbind(rep(1, dim(mu)[1]), mu)) * lambda_hat))
    abs_diff <- dist_euclidean(theta_new, theta_old)                              # The Euclidean distance between the old and new theta vectors.
    if (is.na(abs_diff)) {
      abs_diff <- 0
      theta_new <- theta_old
    } # It helps us avoid "NA trouble".
    
    theta_old <- pmax(theta_new, 0) # Set the new theta as the old theta for the next iteration step.   
    theta_old <- pmin(theta_old, 1) # pmax() and pmin() guarantee that theta_j's are in [0,1]. 
    count <- count + 1
    lambda_hat_old <- lambda_hat
  }
  
  theta_hat <- pmax(theta_new, 0) # The iteration has stopped.  
  theta_hat <- pmin(theta_hat, 1) # Again, we guarantee that theta_j's are in [0,1]
  return(theta_hat) # It returns the estimation of weights \theta_j's_
}

## Subsection 2.2
# Based on the function weight.seq() in Subsection 2.1, we give the following
# high dimensional mixture density estimation function.

hdmde <- function(x_obs, N0, alpha, max_comp) {
  
  # "x_obs" is the data set of interest. 
  # There are n observations, and each observation is a D-dimensional point.
  # x_obs is a n-by-D matrix.
  # "N0" is a predetermined lower bound for N - the number of density components
  # "alpha" is the predetermined confidence level.
  # "max_comp" is the upper bound of the number of components in the desired mixture density.
  
  zalpha <- qnorm(1 - alpha / 2)  
  n_D <- dim(x_obs) # "x.obs" is a n by D matrix. Each row of it denotes a data point in R^D.
  n <- n_D[1]
  D <- n_D[2]
  N <- N0
  
  km <- kmeans(x_obs, N, iter.max = 100, nstart = 100) # Use k-means clustering to get N clusters of x_obs_
  mu <- km$centers # Set the centers of the N clusters as the means of the density components.
  
  sigma_vec <- rep(NA, N) # The following block estimates \sigma_N.
  for(j in 1:N) {                      
    index_temp <- which(km$cluster == j)   
    xi_j <- matrix(x_obs[index_temp, ], nrow = length(index_temp))
    sig_prepare <- function(x) {
      return((dist_euclidean(x, mu[j, ])) ^ 2) 
    }
    s <- apply(xi_j, 1, sig_prepare)
    sigma_vec[j] <- mean(s)
  }
  sig <- sqrt(mean(sigma_vec) / dim(x_obs)[2])
  
  # It gives an estimation of theta_j's with weight.seq().
  theta_hat <- weight_seq(x_obs, mu, sig)  
  
  # The following block gives an approximation to the underlying density function of interest.
  # This estimation is of the form of weights times scaled kernels.
  f_test <- function(x) {
    fun_prepare <- function(t) { 
      return(ker(x, t, sig)) 
    }
    comp_vec <- apply(mu, 1, fun_prepare)
    return(sum(theta_hat * comp_vec))
  }
  
  p_old <- apply(x_obs, 1, f_test) # The second/negative term p_N in Delta.hat.
  
  test_rejection <- 1
  while ((test_rejection == 1) & (N <= min(n, max_comp))) {
    
    N <- N + 1
    
    ##################################################
    # The following is a repetition of the codes above.
    km <- kmeans(x_obs, N, iter.max = 100, nstart = 100)
    mu <- km$centers
    
    sigma_vec <- rep(NA, N)
    for(j in 1:N) {
      index_temp <- which(km$cluster == j)
      xi_j <- matrix(x_obs[index_temp, ], nrow = length(index_temp))
      sig_prepare <- function(x) {
        return(dist_euclidean(x, mu[j, ]) ^ 2)
      }
      s <- apply(xi_j, 1, sig_prepare)
      sigma_vec[j] <- mean(s)
    }
    sig <- sqrt(mean(sigma_vec) / dim(x_obs)[2])
    
    theta_hat <- weight_seq(x_obs, mu, sig) 
    
    f_test <- function(x) {
      fun_prepare <- function(t) {
        return(ker(x, t, sig))
      }
      comp_vec <- apply(mu, 1, fun_prepare)
      return(sum(theta_hat * comp_vec))
    }
    
    # The part above is a repetition.
    ##################################################
    
    p_new <- apply(x_obs, 1, f_test) # The first/positive term p_{N+1} in Delta_hat_
    delta_hat <- p_new - p_old # Delta.hat
    sigma_hat_sq <- mean((delta_hat - mean(delta_hat)) ^ 2)
    Z_I_N <- sqrt(dim(x_obs)[1]) * mean(delta_hat) / sqrt(sigma_hat_sq)
    
    if (
      (Z_I_N <= zalpha) & 
      (Z_I_N >= -zalpha) & 
      (!is.na(Z_I_N))
    ) { 
      test_rejection=0 
    }
    p_old <- p_new
  }
  
  f <- f_test
  
  resp <- list(
    estimating_pdf = f,
    theta_hat = theta_hat,
    mu = mu,
    N = N,
    k_means_result = km,
    sigma = sig,
    Z_I_N = Z_I_N
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

############ Section 3, Principal Manifold Estimation ######################
############################################################################

PME <- function(x_obs, d, N0=20*D, tuning_para_seq=exp((-15:5)), alpha=0.05, max_comp=100, epsilon=0.05, max_iter=100, print_MSDs=TRUE) {
  
  # "x.obs" is the data set of interest. 
  #         There are n observations, and each observation is a D-dimensional point.
  #         x.obs is a n-by-D matrix.
  # "d" is the intrinsic dimension of the underlying manifold
  # "N0" is a predetermined lower bound for N - the number of density components, default value is 20*D
  # "tuning.para.seq" is a vector of tuning parameter candidates, its default value is exp((-15:5)).
  #                   If you would like to fit a manifold for a specific lambda, set tuning.prar.seq=c(lambda).
  # "alpha" is the pre-determined confidence level, which determines the number of the components in a mixture density.
  # "max.comp" is the upper bound of the number of components in the desired mixture density.
  # "epsilon" is the tolerance for distance between manifolds f.old and f.new.
  # "max.iter" is the upper bound of the number of steps in this iteration.
  
  # Remark: The larger N0 is, the less time consuming the function is.
  
  
  dimension_size <- dim(x_obs)
  D <- dimension_size[2] # "D" is the dimension of the input space
  n <- dimension_size[1] # "n" is the number of observed D-dimensional data points
  lambda <- 4 - d # "lambda" determines the form of reproducing kernels
  
  if (N0 == 0) { 
    N0 <- 20 * D 
  }
  
  # 
  est <- hdmde(x_obs, N0, alpha, max_comp) # "hdmde" gives \hat{Q}_N.
  theta_hat <- est$theta_hat
  centers <- est$mu
  sigma <- est$sigma
  W <- diag(theta_hat) # The matrix W
  X <- est$mu
  I <- length(theta_hat)
  
  # The (i,j)th element of this matrix is the Euclidean 
  # distance between mu[i,] and mu[j,].
  dissimilarity_matrix <- as.matrix(dist(X))                               
  # Give the initial projection indices by ISOMAP.
  isomap_initial <- isomap(dissimilarity_matrix, ndim = d, k = 10)            
  t_initial <- isomap_initial$points 
  
  MSE_seq <- vector()
  SOL <- list()
  TNEW <- list()
  
  for (tuning_ind in 1:length(tuning_para_seq)) {
    
    print(
      paste(
        "The tuning parameter is lambda[", 
        as.character(tuning_ind), 
        "] = ", 
        as.character(tuning_para_seq[tuning_ind]), 
        "."
      )
    )
    
    w <- tuning_para_seq[tuning_ind]
    tnew <- t_initial
    t_val <- cbind(rep(1, I), tnew) # The matrix T
    
    E <- matrix(NA, ncol = I, nrow = I)                                          
    for (j in 1:I) {
      temp_mat <- sweep(tnew, 2, tnew[j, ])
      E[, j] <- apply(temp_mat, 1, eta_kernel, lambda)
    }
    
    # This block gives the first step of iteration.
    ###############################################
    M1 <- cbind(
      2 * E %*% W %*% E + 2 * w * E, 
      2 * E %*% W %*% t_val,
      t_val
    )
    M2 <- cbind(
      2 * t(t_val) %*% W %*% E,
      2 * t(t_val) %*% W %*% t_val,
      matrix(0, ncol = d + 1, nrow = d + 1)
    )
    M3 <- cbind(
      t(t_val),
      matrix(0, ncol = d + 1, nrow = d + 1),
      matrix(0, ncol = d + 1, nrow = d + 1)
    )
    M <- rbind(M1, M2, M3) # The coefficient matrix of the linear equations
    
    # The nonhomogeneous term of the linear equations
    b <- rbind(
      2 * E %*% W %*% X,
      2 * t(t_val) %*% W %*% X,
      matrix(0, nrow = d + 1, ncol = D)
    )
    sol <- ginv(M) %*% b # Solve the linear equations
    
    eta_func <- function(t) {
      temp_mat <- sweep(-1 * tnew, 2, t, "+")
      return(
        matrix(
          apply(temp_mat, 1, eta_kernel, lambda),
          ncol = 1
        )
      )
    }

    fnew <- function(t) {
      return(
        as.vector(
          t(sol[1:I, ]) %*% eta_func(t) +
            t(sol[(I + 1):(I + d + 1), ]) %*% matrix(c(1, t), ncol = 1)
        )
      )
    }
    
    # ISOMAP gives the initial projection indices. Then the initial projection indices give the initial manifold f0.
    # The new projection indices are derived by projecting mu_j's onto f0. 
    
    f0 <- fnew
    
    X_initial_guess <- cbind(X, tnew) # The "tnew" here is derived from ISOMAP.
    projection_index_f0 <- function(x_init) { 
      mem_projection(x_init[1:D], f0, x_init[(D + 1):(D + d)]) 
    }
    
    # The first D columns of x.init corresponds to X and the last d columns corresponds to tnew.
    # projection() is applied to X[j,] with initial guess tnew[j,], which is the projection index for X[j,] onto the old manifold f0.
    
    tnew <- matrix(
      t(apply(X_initial_guess, 1, projection_index_f0)),
      nrow = I
    )   # This method can help us avoid the chaos from improper initial guess.
    
    # Sum of the squared distances between x_i and its projection onto manifold f.
    # The first D columns of "x.prin" corresponds to points in the input space
    # and the last d columns of "x.prin" corresponds to the projection indices of these points onto f.
    SSD_prepare <- function(x_prin, f) { 
      return(dist_euclidean(x_prin[1:D], f(x_prin[(D + 1):(D + d)])) ^ 2) 
    }
    
    X_projection_index <- cbind(X, tnew) # "tnew" here is the projection index onto fnew, rather than f0. 
    SSD_prepare_again <- function(x_init) { 
      return(SSD_prepare(x_init, fnew)) 
    }
    
    SSD_new <- apply(X_projection_index, 1, SSD_prepare_again) %>% 
      as.vector() %>% 
      sum()
    
    # This block gives the first step of iteration.
    ###############################################
    
    # The iteration for PME is given by the following loop.
    count <- 1
    SSD_ratio <- 10 * epsilon # A quantity measuring the distance between f0 and fnew.
    
    while (
      (SSD_ratio > epsilon) & 
      (SSD_ratio <= 5) & 
      (count <= (max_iter - 1))
    ) {
      
      SSD_old <- SSD_new
      f0 <- fnew # Set the new manifold in the previous step as the old manifold in the next step.
      
      # Repetition 
      #############
      
      t_val <- cbind(rep(1, I), tnew)                                             
      
      E <- matrix(NA, ncol = I, nrow = I)                                          
      for (j in 1:I) {
        temp_mat <- sweep(tnew, 2, tnew[j, ])
        E[, j] <- apply(temp_mat, 1, eta_kernel, lambda)
      }
      
      # This block gives the first step of iteration.
      ###############################################
      M1 <- cbind(
        2 * E %*% W %*% E + 2 * w * E,
        2 * E %*% W %*% t_val,
        t_val
      )
      M2 <- cbind(
        2 * t(t_val) %*% W %*% E,
        2 * t(t_val) %*% W %*% t_val, 
        matrix(0, ncol = d + 1, nrow = d + 1)
      )
      M3 <- cbind(
        t(t_val),
        matrix(0, ncol = d + 1, nrow = d + 1),
        matrix(0, ncol = d + 1, nrow = d + 1)
      )
      M <- rbind(M1, M2, M3)
      
      b <- rbind(
        2 * E %*% W %*% X,
        2 * t(t_val) %*% W %*% X,
        matrix(0, nrow = d + 1, ncol = D))      
      sol <- ginv(M) %*% b
      
      eta_func <- function(t) {
        temp_mat <- sweep(-1 * tnew, 2, t, "+")
        return(
          matrix(
            apply(temp_mat, 1, eta_kernel, lambda),
            ncol = 1
          )
        )
      }
      fnew <- function(t) {
        return(
          as.vector(
            t(sol[1:I, ]) %*% eta_func(t) +
              t(sol[(I + 1):(I + d + 1), ]) %*% matrix(c(1, t), ncol = 1)
          )
        )
      }
      
      t_old <- tnew
      # Repetition 
      #############
      
      X_initial_guess <- cbind(X, tnew) # The "tnew" here is the projection index to f0.
      projection_index_f0 <- function(x_init) { 
        mem_projection(x_init[1:D], f0, x_init[(D + 1):(D + d)]) 
      }
      
      tnew <- matrix(
        t(apply(X_initial_guess, 1, projection_index_f0)),
        nrow = I
      )
      
      X_projection_index <- cbind(X, tnew)
      SSD_prepare_again <- function(x_init) { 
        return(SSD_prepare(x_init, fnew)) 
      }
      SSD_new <- apply(
        X_projection_index,
        1,
        SSD_prepare_again
      ) %>% 
        as.vector() %>% 
        sum()
      
      SSD_ratio <- abs(SSD_new - SSD_old) / SSD_old
      count <- count + 1
      
      print(
        paste(
          "SSD.ratio is ",
          as.character(round(SSD_ratio, 4)),
          " and this is the ",
          as.character(count),
          "th step of iteration."
        )
      )
    }
    
    
    # For a fixed tuning parameter value, the corresponding MSD is computed by the following chunk.
    km <- est$k_means_result
    data_initial <- matrix(0, nrow = 1, ncol = D + d)
    for(i in 1:I) {
      index_temp <- which(km$cluster == i)
      length_temp <- length(index_temp)
      X_i <- matrix(x_obs[index_temp, ], nrow = length_temp)
      t_temp <- matrix(rep(tnew[i, 1], length_temp))
      for(j in 1:d) {
        t_temp <- cbind(t_temp, rep(tnew[i, j], length_temp))
      }
      t_temp <- matrix(t_temp[, -1], nrow = length_temp)
      data_initial <- rbind(data_initial, cbind(X_i, t_temp))
    }
    data_initial <- data_initial[-1, ]
    proj_para_prepare <- function(data_init) { 
      return(projection(data_init[1:D], fnew, data_init[(D+1):(D+d)])) 
    }
    proj_para <- matrix(
      t(apply(data_initial, 1, proj_para_prepare)),
      ncol = d
    )
    proj_points <- t(apply(proj_para, 1, fnew))
    diff_data_fit <- apply(
      data_initial[, 1:D] - proj_points,
      1,
      norm_euclidean
    )
    MSE <- mean(diff_data_fit ^ 2)
    
    MSE_seq[tuning_ind] <- MSE
    print(
      paste(
        "When lambda = ",
        as.character(w),
        ", MSD = ",
        as.character(MSE),
        "."
      )
    )
    SOL[[tuning_ind]] <- sol
    TNEW[[tuning_ind]] <- tnew
    
    # To reduce the computational burden, if the MSD in the k-th step of the for-loop is
    # smaller than that in the next 4 steps of this for-loop (k+1, k+2, k+3, k+4), 
    # we stop this for-loop. 
    if (tuning_ind >= 4) {
      if (
        (MSE_seq[tuning_ind] > MSE_seq[tuning_ind - 1]) & 
        (MSE_seq[tuning_ind - 1] > MSE_seq[tuning_ind - 2]) & 
        (MSE_seq[tuning_ind - 2] > MSE_seq[tuning_ind - 3])
      ) {
        break
      }
    }
  }
  
  # The following chunk gives the f_\lambda with the optimal \lambda.
  optimal_ind <- min(which(MSE_seq == min(MSE_seq)))
  sol_opt <- SOL[[optimal_ind]]
  tnew_opt <- TNEW[[optimal_ind]]
  eta_func <- function(t) {
    eta_func_prepare <- function(tau) { 
      return(eta_kernel(t - tau, lambda)) 
    }
    return(matrix(apply(tnew_opt, 1, eta_func_prepare), ncol = 1))
  }
  
  f_optimal <- function(t) {
    return(
      as.vector(
        t(sol_opt[1:I, ]) %*%
          eta_func(t) + t(sol_opt[(I + 1):(I + d + 1), ]) %*%
          matrix(c(1, t), ncol=1)
      )
    )
  }
  
  
  if (print_MSDs == TRUE) {
    plot(
      log(tuning_para_seq[1:length(MSE_seq)]), 
      MSE_seq, 
      xlab = "Log Lambda", 
      ylab = "MSD", 
      type = "l"
    )
    lines(
      log(tuning_para_seq[1:length(MSE_seq)]), 
      MSE_seq, 
      type = "p",
      pch = 20,
      col = "orange",
      cex = 2
    )
    abline(
      v = log(tuning_para_seq[optimal_ind]),
      lwd = 1.5,
      col = "darkgreen",
      lty = 2
    )
    
    print(
      paste(
        "The optimal tuning parameter is ", 
        as.character(tuning_para_seq[optimal_ind]), 
        ", and the MSD of the optimal fit is ",
        as.character(MSE_seq[optimal_ind]),
        "."
      )
    )
  }
  
  resp <- list(
    embedding_map = f_optimal, 
    MSD = MSE_seq,  
    knots = centers,
    weights_of_knots = theta_hat,
    coe_kernel = sol_opt[1:I, ],
    coe_poly = sol_opt[(I + 1):(I + d + 1), ],
    SOL = SOL,
    TNEW = TNEW,
    T_parameter = sol_opt,
    Coef = tnew_opt
  )
  return(resp)
  
  # Output:
  # embedding.map: The fitted embedding map R^d -> R^D. 
  # MSD: A vector of mean squared distances. 
  #      Each component of this vector corresponds to a tuning parameter candidate.
  # knots: Knots in the discrete measure \hat{Q}_N
  # weights.of.knots: A vector of weights for the knots of \hat{Q}_N. 
  #                   The k-th component is the weight for the k-th knot.
  # Lists T.parameter and Coef are quantities determining the analytic formula of f.optimal.
  # coe.kernel and coe.poly are quantities for the Interior Identification function.
  
}

##################################################################
## Section 3 completes the principal manifold estiamtion function.
##################################################################
