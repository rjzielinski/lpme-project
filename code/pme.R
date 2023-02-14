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
library("MASS")
#install.packages("Matrix")
library("Matrix")
#install.packages("vegan")
library("vegan")
#install.packages("plot3D")
library("plot3D")
library(Rfast)
library(tidyverse)
library(plotly)
library(Rcpp)


########### Section 1, Some Basic Functions ########################
####################################################################

## Subsection 1.1, Functions for Euclidean metrics
source("code/functions/norm_euclidean.R")

# Distance function in an Euclidean space of any dimension
source("code/functions/dist_euclidean.R")

## Subsection 1.2, Kernels for minimization in a semi-normed space of Sobolev type

source("code/functions/ker.R")

source("code/functions/eta_kernel.R")

## Subsection 1.3, Projection Index function

source("code/functions/projection.R")

##### Section 2, High Dimensional Mixture Density Estimation #######
####################################################################

## Subsection 2.1
# When \mu's and \sigma are given, the following function estimates \hat{\theta}'s.

source("code/functions/weight_seq.R")

## Subsection 2.2
# Based on the function weight.seq() in Subsection 2.1, we give the following
# high dimensional mixture density estiamtion function.

source("code/functions/hdmde.R")

source("code/functions/solve_eq.R")

sourceCpp("code/functions/pme_functions.cpp")

############ Section 3, Principal Manifold Estimation ######################
############################################################################

pme <- function(x.obs, d, initialization = NULL, N0=20*D, tuning.para.seq=exp((-15:5)), alpha=0.05, max.comp=100, epsilon=0.05, max.iter=100, print.MSDs=TRUE, print_plots = FALSE) {

  # "x.obs" is the data set of interest.
  #         There are n observations, and each observation is a D-dimensional point.
  #         x.obs is a n-by-D matrix.
  # "d" is the intrinsic dimension of the underlying manifold
  # "initialization" is an optional parameter in the form of a list with
  # an output from the hdmde function in the first index and isomap object in the second
  # "N0" is a predetermined lower bound for N - the number of density components, default value is 20*D
  # "tuning.para.seq" is a vector of tuning parameter candidates, its default value is exp((-15:5)).
  #                   If you would like to fit a manifold for a specific lambda, set tuning.prar.seq=c(lambda).
  # "alpha" is the pre-determined confidence level, which determines the number of the components in a mixture density.
  # "max.comp" is the upper bound of the number of components in the desired mixture density.
  # "epsilon" is the tolerance for distance between manifolds f.old and f.new.
  # "max.iter" is the upper bound of the number of steps in this iteration.

  # Remark: The larger N0 is, the less time consuming the function is.


  dimension.size <- dim(x.obs)
  D <- dimension.size[2] # "D" is the dimension of the input space.
  n <- dimension.size[1] # "n" is the number of observed D-dimensional data points.
  lambda <- 4 - d # "lambda" determines the form of reproducing kernels

  if (N0 == 0) {
    N0 <- 20 * D
  }

  if (is.null(initialization)) {
    est <- hdmde(x.obs, N0, alpha, max.comp) # "hdmde" gives \hat{Q}_N.
    est_order <- order(est$mu[, 1])
    theta.hat <- est$theta.hat[est_order]
    centers <- est$mu[est_order, ]
    sigma <- est$sigma
    W <- diag(theta.hat) # The matrix W
    X <- est$mu[est_order, ]
    I <- length(theta.hat)

    # The (i,j)th element of this matrix is the Euclidean
    # distance between mu[i,] and mu[j,].
    dissimilarity.matrix <- as.matrix(dist(X))
    isomap.initial <- isomap(dissimilarity.matrix, ndim = d, k = 10)
    t.initial <- isomap.initial$points
  } else {
    isomap.initial <- initialization[[1]]
    theta.hat <- initialization[[2]]
    centers <- initialization[[3]]
    sigma <- initialization[[4]]
    W <- diag(theta.hat)
    X <- centers
    I <- length(theta.hat)
    t.initial <- initialization[[5]]
  }

  MSE.seq <- vector()
  SOL <- list()
  TNEW <- list()

  for (tuning.ind in 1:length(tuning.para.seq)) {

    print(
      paste(
        "The tuning parameter is lambda[",
        as.character(tuning.ind),
        "] = ",
        as.character(tuning.para.seq[tuning.ind]),
        "."
      )
    )

    w <- tuning.para.seq[tuning.ind]
    tnew <- t.initial
    t_val <- cbind(rep(1, I), tnew) # The matrix T

    # E <- matrix(NA, ncol = I, nrow = I)
    # for(j in 1:I) {
    #   E.prepare <- function(t) {
    #     eta_kernel(t - tnew[j, ], lambda)
    #   }
    #   E[, j] <- apply(tnew, 1, E.prepare) # The matrix E
    # }

    E <- calcE(tnew, lambda)

    # This block gives the first step of iteration.
    ###############################################


    sol <- solve_eq(E, W, t_val, X, w, d, D)

    # eta.func <- function(t) {
    #   eta.func.prepare <- function(tau) {
    #     return(eta_kernel(t - tau, lambda))
    #   }
    #   return(
    #     matrix(
    #       apply(tnew, 1, eta.func.prepare),
    #       ncol=1
    #     )
    #   )
    # }

    fnew <- function(t) {
      return(
        as.vector(
          (t(sol[1:I, ]) %*% etaFunc(t, tnew, lambda)) +
            (t(sol[(I + 1):(I + d + 1), ]) %*% matrix(c(1, t), ncol = 1))
        )
      )
    }

    fnew2 <- function(t, sol, tnew, I, d, lambda) {
      out_mat <- (t(sol[1:I, ]) %*% etaFunc2(t, tnew, lambda)) +
          (t(sol[(I + 1):(I + d + 1), ]) %*% t(cbind(1, t)))
      return(t(out_mat))
    }

    # ISOMAP gives the initial projection indices. Then the initial projection indices give the initial manifold f0.
    # The new projection indices are derived by projecting mu_j's onto f0.

    f0 <- fnew

    X.initial.guess <- cbind(X, tnew) # The "tnew" here is derived from ISOMAP.
    projection.index.f0 <- function(x.init) {
      projection(x.init[1:D], f0, x.init[(D + 1):(D + d)])
    }

    # The first D columns of x.init corresponds to X and the last d columns corresponds to tnew.
    # projection() is applied to X[j,] with initial guess tnew[j,], which is the projection index for X[j,] onto the old manifold f0.

    # tnew <- map(
    #   1:I,
    #   ~ projection(
    #     X.initial.guess[.x, 1:D],
    #     f0,
    #     X.initial.guess[.x, (D + 1):(D + d)]
    #   )
    # ) %>%
    #   unlist() %>%
    #   matrix(nrow = I, byrow = TRUE)

    tnew <- calc_tnew(X, tnew, sol, I, d, lambda)

    # Sum of the squared distances between x_i and its projection onto manifold f.
    # The first D columns of "x.prin" corresponds to points in the input space
    # and the last d columns of "x.prin" corresponds to the projection indices of these points onto f.
    X.projection.index <- cbind(X, tnew)
    SSD.new <- map(
      1:I,
      ~ dist_euclideanC(
        X.projection.index[.x, 1:D],
        fNew(X.projection.index[.x, (D + 1):(D + d)], sol, tnew, I, d, lambda)
      ) ^ 2
    ) %>%
      unlist() %>%
      sum()

    # This block gives the first step of iteration.
    ###############################################
    if (print_plots == TRUE) {
      plot_pme(x.obs, sol, tnew, I, d, lambda)
    }

    # The iteration for PME is given by the following loop.
    count <- 1
    SSD.ratio <- 10 * epsilon # A quantity measuring the distance between f0 and fnew.

    while (
      (SSD.ratio > epsilon) &
      (SSD.ratio <= 5) &
      (count <= (max.iter - 1))
    ) {

      SSD.old <- SSD.new
      f0 <- fnew # Set the new manifold in the previous step as the old manifold in the next step.

      # Repetition
      #############

      t_val <- cbind(rep(1, I), tnew)

      E <- calcE(tnew, lambda)

      sol <- solve_eq(E, W, t_val, X, w, d, D)

      fnew <- function(t) {
        as.vector(
          t(sol[1:I, ]) %*%
            etaFunc(t, tnew, lambda) + t(sol[(I + 1):(I + d + 1), ]) %*%
            # eta.func(t) + t(sol[(I + 1):(I + d + 1), ]) %*%
            matrix(c(1, t), ncol = 1)
        ) %>%
          return()
      }


      t.old <- tnew
      # Repetition
      #############

      # X.initial.guess <- cbind(X, tnew) # The "tnew" here is the projection index to f0.
      # projection.index.f0 <- function(x.init) {
      #   projection(x.init[1:D], f0, x.init[(D + 1):(D + d)])
      # }
      # tnew <- matrix(
      #   t(apply(X.initial.guess, 1, projection.index.f0)),
      #   nrow = I
      # ) # The "tnew" here is the projection index to fnew, rather than f0.

      tnew <- calc_tnew(X, tnew, sol, I, d, lambda)

      X.projection.index <- cbind(X, tnew)
      # SSD.prepare.again <- function(x.init) {
      #   return(SSD.prepare(x.init, fnew))
      # }
      # SSD.new <- apply(
      #   X.projection.index,
      #   1,
      #   SSD.prepare.again
      # ) %>%
      #   as.vector() %>%
      #   sum()

      SSD.new <- map(
        1:I,
        ~ dist_euclideanC(
          X[.x, ],
          fNew(tnew[.x, ], sol, tnew, I, d, lambda)
        ) ^ 2
      ) %>%
        unlist() %>%
        sum()


      if (print_plots == TRUE) {
        plot_pme(x.obs, sol, tnew, I, d, lambda)
      }

      SSD.ratio <- abs(SSD.new - SSD.old) / SSD.old
      count <- count + 1

      print(
        paste0(
          "SSD is ",
          as.character(round(SSD.new, 4)),
          " SSD.ratio is ",
          as.character(round(SSD.ratio, 4)),
          " and this is the ",
          as.character(count),
          "th step of iteration."
        )
      )
    }

    # For a fixed tuning parameter value, the corresponding MSD is computed by the following chunk.
    if (is.null(initialization)) {
      km <- est$k.means.result
    } else {
      km <- initialization[[6]]
    }
    data.initial <- matrix(0, nrow = 1, ncol = D + d)
    for(i in 1:I) {
      index.temp <- which(km$cluster == i)
      length.temp <- length(index.temp)
      X.i <- matrix(x.obs[index.temp, ], nrow = length.temp)
      t.temp <- matrix(rep(tnew[i, 1], length.temp))
      for(j in 1:d) {
        t.temp <- cbind(t.temp, rep(tnew[i, j], length.temp))
      }
      t.temp <- matrix(t.temp[, -1], nrow = length.temp)
      data.initial <- rbind(data.initial, cbind(X.i, t.temp))
    }
    data.initial <- data.initial[-1, ]
    proj.para.prepare <- function(data.init) {
      return(projection(data.init[1:D],fnew,data.init[(D+1):(D+d)]))
    }
    proj.para <- matrix(
      t(apply(data.initial, 1, proj.para.prepare)),
      ncol = d
    )
    proj.points <- t(apply(proj.para, 1, fnew))
    diff.data.fit <- apply(
      data.initial[, 1:D] - proj.points,
      1,
      norm_euclidean
    )
    MSE <- mean(diff.data.fit ^ 2)

    MSE.seq[tuning.ind] <- MSE
    print(
      paste(
        "When lambda = ",
        as.character(w),
        ", MSD = ",
        as.character(MSE),
        "."
      )
    )
    SOL[[tuning.ind]] <- sol
    TNEW[[tuning.ind]] <- tnew

    # To reduce the computational burden, if the MSD in the k-th step of the for-loop is
    # smaller than that in the next 4 steps of this for-loop (k+1, k+2, k+3, k+4),
    # we stop this for-loop.
    if (tuning.ind >= 4) {
      if (
        (MSE.seq[tuning.ind] > MSE.seq[tuning.ind - 1]) &
        (MSE.seq[tuning.ind - 1] > MSE.seq[tuning.ind - 2]) &
        (MSE.seq[tuning.ind - 2] > MSE.seq[tuning.ind - 3])
      ) {
        break
      }
    }
  }

  # The following chunk gives the f_\lambda with the optimal \lambda.
  optimal.ind <- min(which(MSE.seq == min(MSE.seq)))
  sol.opt <- SOL[[optimal.ind]]
  tnew.opt <- TNEW[[optimal.ind]]
  # eta.func <- function(t) {
  #   eta.func.prepare <- function(tau) {
  #     return(eta_kernel(t - tau, lambda))
  #   }
  #   return(matrix(apply(tnew.opt, 1, eta.func.prepare), ncol = 1))
  # }

  f.optimal <- function(t) {
    return(
      as.vector(
        t(sol.opt[1:I, ]) %*%
          etaFunc(t, tnew.opt, lambda) + t(sol.opt[(I + 1):(I + d + 1), ]) %*%
          # eta.func(t) + t(sol.opt[(I + 1):(I + d + 1), ]) %*%
          matrix(c(1, t), ncol=1)
      )
    )
  }


  if (print.MSDs == TRUE) {
    plot(
      log(tuning.para.seq[1:length(MSE.seq)]),
      MSE.seq,
      xlab = "Log Lambda",
      ylab = "MSD",
      type = "l"
    )
    lines(
      log(tuning.para.seq[1:length(MSE.seq)]),
      MSE.seq,
      type = "p",
      pch = 20,
      col = "orange",
      cex = 2
    )
    abline(
      v = log(tuning.para.seq[optimal.ind]),
      lwd = 1.5,
      col = "darkgreen",
      lty = 2
    )

    print(
      paste(
        "The optimal tuning parameter is ",
        as.character(tuning.para.seq[optimal.ind]),
        ", and the MSD of the optimal fit is ",
        as.character(MSE.seq[optimal.ind]),
        "."
      )
    )
  }

  resp <- list(
    embedding.map = f.optimal,
    MSD = MSE.seq,
    knots = km,
    weights.of.knots = theta.hat,
    coe.kernel = sol.opt[1:I, ],
    coe.poly = sol.opt[(I + 1):(I + d + 1), ],
    SOL = SOL,
    TNEW = TNEW,
    T.parameter = sol.opt,
    Coef = tnew.opt
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
