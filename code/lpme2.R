lpme <- function(df, d, tuning.para.seq = exp(-15:5), alpha = 0.05, max.comp = 500, epsilon = 0.05, max.iter = 100, verbose = "all", print_plots = TRUE, SSD_ratio_threshold = 100, increase_threshold = 1.05, init = "full") {
  # df is an N x (D + 1) matrix, with the first column corresponding
  # to the time point at which each observation was collected
  # this matrix should include the observations from all time points

  source("code/pme.R")
  source("code/functions/plot_lpme.R")
  source("code/functions/solve_eq_inter_t.R")
  source("code/functions/projection_lpme.R")
  source("code/functions/solve_eq2.R")
  require(plotly)
  require(svMisc)
  require(Rcpp)
  sourceCpp("code/functions/pme_functions.cpp")
  # source("Principal_Manifold_Estimation.R")

  if (verbose == "none") {
    print.MSDs <- FALSE
    print.SSDs <- FALSE
  } else if (verbose == "MSD") {
    print.MSDs <- TRUE
    print.SSDs <- FALSE
  } else if (verbose == "SSD") {
    print.MSDs <- FALSE
    print.SSDs <- TRUE
  } else if (verbose == "all") {
    print.MSDs <- TRUE
    print.SSDs <- TRUE
  }

  time_points <- df[, 1] %>%
    unique()

  pme_results <- list()
  pme_coefs <- list()
  pme_coefs2 <- list()
  funcs <- list()
  clusters <- list()
  centers <- list()
  x_test <- list()
  r <- list()
  r2 <- list()

  init_timevals <- list()
  init_theta_hat <- list()
  init_centers <- list()
  init_sigma <- list()
  init_clusters <- list()


  if (init %in% c("first", "full")) {
    if (init == "first") {
      init_df <- df[df[, 1] == time_points[1], -1]
    } else if (init == "full") {
      for (idx in 1:length(time_points)) {
        init_dimension_size <- dim(df[, -1])
        init_D <- init_dimension_size[2]
        init_n <- init_dimension_size[1]
        lambda <- 4 - d
        init_N0 <- 20 * init_D
        init_df_temp <- df[df[, 1] == time_points[idx], -1]
        init_est_temp <- hdmde(init_df_temp, init_N0, alpha, max.comp)
        init_timevals[[idx]] <- rep(time_points[idx], dim(init_est_temp$mu)[1])
        init_theta_hat[[idx]] <- init_est_temp$theta.hat
        init_centers[[idx]] <- cbind(time_points[idx], init_est_temp$mu)
        init_sigma[[idx]] <- init_est_temp$sigma
        init_clusters[[idx]] <- init_est_temp$k.means.result
      }
    }

    init_timevals <- reduce(init_timevals, c)
    init_theta_hat <- reduce(init_theta_hat, c)
    init_centers <- reduce(init_centers, rbind)
    init_sigma <- reduce(init_sigma, c)
    init_W <- diag(init_theta_hat)
    init_X <- init_centers[, -1]
    init_I <- length(init_theta_hat)

    init_dissimilarity_matrix <- as.matrix(dist(init_X))
    init_isomap <- isomap(init_dissimilarity_matrix, ndim = d, k = 10)

    # init_dimension.size <- dim(init_df)
    # init_D <- init_dimension.size[2]
    # init_n <- init_dimension.size[1]
    # lambda <- 4 - d
    #
    # init_N0 <- 20 * init_D
    # init_est <- hdmde(init_df, init_N0, alpha, max.comp)
    # init_order <- order(init_est$mu[, 1])
    # init_theta.hat <- init_est$theta.hat[init_order]
    # init_centers <- init_est$mu[init_order, ]
    # init_sigma <- init_est$sigma
    # init_W <- diag(init_theta.hat)
    # init_X <- init_est$mu[init_order, ]
    # init_I <- length(init_theta.hat)
    #
    # init_dissimilarity.matrix <- as.matrix(dist(init_X))
    # init_isomap <- isomap(init_dissimilarity.matrix, ndim = d, k = 10) # could changing k improve fit?
  }

  num_clusters <- rep(0, length(time_points))

  for (idx in 1:length(time_points)) {
    df_temp <- df[df[, 1] == time_points[idx], ]
    if (init %in% c("first", "full")) {
      pme_results[[idx]] <- pme(
        # x.obs = df_temp[, -1],
        x.obs = df_temp,
        d = d + 1,
        initialization = list(
          init_isomap,
          init_theta_hat[init_timevals == time_points[idx]],
          init_centers[init_timevals == time_points[idx], ],
          init_sigma[idx],
          cbind(
            init_centers[init_timevals == time_points[idx], 1],
            matrix(init_isomap$points[init_timevals == time_points[idx], ], nrow = length(init_theta_hat[init_timevals == time_points[idx]]))
          ),
          init_clusters[[idx]]
        ),
        verbose = "MSD"
      )
    } else {
      pme_results[[idx]] <- pme(
        x.obs = df_temp[, -1],
        d = d,
        verbose = "MSD"
      )
    }
    opt_run <- which.min(pme_results[[idx]]$MSD)
    funcs[[idx]] <- pme_results[[idx]]$embedding.map
    if (idx == 1) {
      centers[[idx]] <- pme_results[[idx]]$knots$centers
      clusters[[idx]] <- pme_results[[idx]]$knots$cluster
      num_clusters[idx] <- dim(pme_results[[idx]]$knots$centers)[1]
    } else {
      centers[[idx]] <- pme_results[[idx]]$knots$centers
      clusters[[idx]] <- pme_results[[idx]]$knots$cluster + sum(num_clusters)
      num_clusters[idx] <- dim(pme_results[[idx]]$knots$centers)[1]
    }

    pme_coefs[[idx]] <- pme_results[[idx]]$SOL[[opt_run]][1:num_clusters[idx], ]
    pme_coefs2[[idx]] <- pme_results[[idx]]$SOL[[opt_run]][(num_clusters[idx] + 1):(num_clusters[idx] + d + 2), ] %>%
      t() %>%
      matrix(nrow = 1)
    # pme_coefs2[[idx]] <- pme_results[[idx]]$SOL[[opt_run]][(num_clusters[idx] + 1):(num_clusters[idx] + d + 2), ]


    r_mat <- pme_results[[idx]]$TNEW[[opt_run]]
    # x_test[[idx]] <- apply(
    #   r_mat,
    #   1,
    #   funcs[[idx]]
    # ) %>%
    #   matrix(nrow = dim(r_mat)[1], byrow = TRUE)
    x_temp <- apply(
      r_mat,
      1,
      funcs[[idx]]
    ) %>%
      t()
    x_test[[idx]] <- x_temp
    # x_test[[idx]] <- x_temp

    r[[idx]] <- cbind(time_points[idx], r_mat)[, -2]
    r2[[idx]] <- time_points[idx]
  }

  r_full <- reduce(r, rbind)
  r_full2 <- reduce(r2, rbind)
  r_df <- data.frame(r_full)
  r_df2 <- data.frame(r_full2)
  rnames <- paste0("r", 1:(dim(r_full)[2] - 1))
  names(r_df) <- c("time", rnames)
  x_test <- reduce(x_test, rbind)
  x_test_df <- data.frame(x_test)
  xnames <- paste0("x", 1:(dim(x_test)[2]))
  names(x_test_df) <- xnames

  coef_full <- reduce(pme_coefs, rbind)
  coef_full2 <- reduce(pme_coefs2, rbind)

  centers_full <- reduce(centers, rbind)
  centers_full <- cbind(x_test[, 1], centers_full)
  clusters_full <- unlist(clusters)

  D_coef <- dim(coef_full)[2]
  D_coef2 <- dim(coef_full2)[2]
  D_new <- dim(x_test)[2]
  d_new <- dim(r_full)[2]
  d_new2 <- dim(r_full2)[2]
  n <- dim(x_test)[1]
  gamma <- 4 - d_new
  gamma2 <- 4 - d_new2

  sig_new <- 0.001
  # theta_hat_new <- weight_seq(x_test, x_test, sig_new)
  theta_hat_new <- init_theta_hat / sum(init_theta_hat)
  centers_new <- x_test
  sigma_new <- 0.001
  W_new <- diag(theta_hat_new)
  W_new2 <- diag(rep(1 / nrow(r_full2), nrow(r_full2)))
  X_new <- x_test
  I_new <- length(theta_hat_new)

  dissimilarity_matrix_new <- as.matrix(dist(X_new))
  t_initial <- r_full %>%
    as.matrix()

  MSE_seq_new <- vector()
  SOL_new <- list()
  TNEW_new <- list()

  for (tuning.ind in 1:length(tuning.para.seq)) {
    if (verbose == "all") {
      print(
        paste0(
          "The tuning parameter is gamma[",
          as.character(tuning.ind),
          "] = ",
          as.character(tuning.para.seq[tuning.ind]),
          "."
        )
      )
    }


    w_new <- tuning.para.seq[tuning.ind]
    t_new <- t_initial
    T_new <- cbind(rep(1, I_new), t_new)
    T_new2 <- cbind(rep(1, nrow(r_full2)), r_full2)

    # E_new <- matrix(NA, ncol = I_new, nrow = I_new)
    # for (j in 1:I_new) {
    #   E_prepare <- function(t) {
    #     eta_kernel(t - t_new[j, ], gamma)
    #   }
    #   E_new[, j] <- apply(t_new, 1, E_prepare)
    # }

    E_new <- calcE(t_new, gamma)
    E_new2 <- calcE(r_full2, gamma2)

    # M1_new <- cbind(
    #   2 * E_new %*% W_new %*% E_new + 2 * w_new * E_new,
    #   2 * E_new %*% W_new %*% T_new,
    #   T_new
    # )
    # M2_new <- cbind(
    #   2 * t(T_new) %*% W_new %*% E_new,
    #   2 * t(T_new) %*% W_new %*% T_new,
    #   matrix(0, ncol = d_new + 1, nrow = d_new + 1)
    # )
    # M3_new <- cbind(
    #   t(T_new),
    #   matrix(0, ncol = d_new + 1, nrow = d_new + 1),
    #   matrix(0, ncol = d_new + 1, nrow = d_new + 1)
    # )
    # M_new <- rbind(
    #   M1_new,
    #   M2_new,
    #   M3_new
    # )
    #
    # b_new <- rbind(
    #   2 * E_new %*% W_new %*% X_new,
    #   2 * t(T_new) %*% W_new %*% X_new,
    #   matrix(0, nrow = d_new + 1, ncol = D_new)
    # )
    # sol_new <- ginv(M_new) %*% b_new

    # sol_new <- solve_eq_inter_t(E_new, W_new, T_new, X_new, w_new, d_new, D_new)
    # sol_new <- solve_eq(E_new, W_new, T_new, X_new, w_new, d_new, D_new)
    # sol_coef <- solve_eq(E_new, W_new, T_new, coef_full, w_new, d_new, D_coef)
    # sol_coef2 <- solve_eq(E_new2, W_new2, T_new2, coef_full2, w_new, d_new2, D_coef2)
    sol_coef <- solve_eq2(E_new, T_new, coef_full, I_new, w_new, d_new, D_coef)
    # sol_coef2 <- solve_eq2(E_new2, T_new2, coef_full2, nrow(r_full2), w_new, d_new2, D_coef2)
    T_new_idx <- lapply(T_new[, 2], function(x) which(time_points == x)) %>%
      reduce(c)
    sol_coef2 <- solve_eq2(E_new, T_new, coef_full2[T_new_idx, ], I_new, w_new, d_new, D_coef2)

    # eta.func <- function(t) {
    #   eta.func.prepare <- function(tau) {
    #     return(eta_kernel(t - tau, gamma))
    #   }
    #   return(
    #     matrix(
    #       apply(
    #         t_new,
    #         1,
    #         eta.func.prepare
    #       ),
    #       ncol = 1
    #     )
    #   )
    # }

    f_coef <- function(t) {
      return_vec <- as.vector(
        t(sol_coef[1:I_new, ]) %*% etaFunc(t, t_new, gamma) +
          t(sol_coef[(I_new + 1):(I_new + d_new + 1), ]) %*% matrix(c(1, t), ncol = 1)
      )
      return(return_vec)
    }

    # f_coef2 <- function(t) {
    #   return_vec <- as.vector(
    #     t(sol_coef2[1:nrow(r_full2), ]) %*% etaFunc(t, r_full2, gamma2) +
    #       t(sol_coef2[(nrow(r_full2) + 1):(nrow(r_full2) + d_new2 + 1), ]) %*% matrix(c(1, t), ncol = 1)
    #   )
    #   return_mat <- matrix(return_vec, ncol = D_coef, byrow = TRUE)
    #   return(return_mat)
    # }
    f_coef2 <- function(t) {
      return_vec <- as.vector(
        t(sol_coef2[1:I_new, ]) %*% etaFunc(t, t_new, gamma) +
          t(sol_coef2[(I_new + 1):(I_new + d_new + 1), ]) %*% matrix(c(1, t), ncol = 1)
      )
      return(return_vec)
    }

    sol_new <- apply(t_new, 1, f_coef) %>%
      t()
    sol_new2 <- apply(t_new, 1, f_coef2) %>%
      t()

    f_new <- function(t) {
      return_vec <- as.vector(
        t(sol_new) %*% etaFunc(t, t_new, gamma) +
          # etaFunc(t, t_new, gamma) + t(sol_new2[(I_new + 1):(I_new + d_new + 1),]) %*%
          # t(f_coef2(t[1])) %*% matrix(c(1, t[1]), ncol = 1)
          t(matrix(f_coef2(t), nrow = d_new + 1, byrow = TRUE)) %*% matrix(c(1, t), ncol = 1)
        )
      # return_vec[1] <- t[1]
      return(return_vec)
    }
    # f_new <- function(t) {
    #   return_vec <- as.vector(
    #     t(coef_full) %*% etaFunc(t, t_new, gamma) +
    #       t(matrix(coef_full2[T_new_idx, ], ncol = D_coef, byrow = TRUE)) %*% matrix(c(1, t), ncol = 1)
    #   )
    #   return_vec
    # }

    # f_new <- function(t) {
    #   return_vec <- as.vector(
    #     t(sol_new[1:I_new, ]) %*%
    #       etaFunc(t, t_new, gamma) + t(sol_new[(I_new + 1):(I_new + d_new + 1),]) %*%
    #       matrix(c(1, t), ncol = 1)
    #     )
    #   # return_vec[1] <- t[1]
    #   return(return_vec)
    # }

    f0_new <- f_new

    X_initial_guess <- cbind(X_new, t_new)
    projection.index.f0 <- function(x.init) {
      projection_lpme(
        x.init[1:D_new],
        f0_new,
        x.init[(D_new + 1):(D_new + d_new)]
      )
    }

    t_new <- matrix(
      t(apply(X_initial_guess, 1, projection.index.f0)),
      nrow = I_new
    )

    # t_new <- calc_tnew(X_new, t_new, sol_new, I_new, d_new, gamma)
    # t_new[, 1] <- t_initial[, 1]

  # SSD.prepare <- function(x.prin, f) {
  #   return(
  #     dist_euclidean(
  #       x.prin[1:D_new],
  #       f(x.prin[(D_new + 1):(D_new + d_new)])
  #     ) ^ 2
  #   )
  # }

    X_projection_index <- cbind(X_new, t_new)
    # SSD.prepare.again <- function(x.init) {
    #   return(SSD.prepare(x.init, f_new))
    # }
    # SSD_new <- apply(X_projection_index, 1, SSD.prepare.again) %>%
    #   as.vector() %>%
    #   sum()

    SSD_new <- map(
      1:I_new,
      ~ dist_euclideanC(
        X_new[.x, ],
        # fNew(t_new[.x, ], sol_new, t_new, I_new, d_new, gamma)
        f_new(t_new[.x, ])
      ) ^ 2
    ) %>%
      unlist() %>%
      sum()

    if (print_plots == TRUE) {
      plot_lpme(df, f_new, t_new, d_new, D_new, time_points)
    }

    count <- 1
    SSD_ratio <- 10 * epsilon

    while ((SSD_ratio > epsilon) & (SSD_ratio <= SSD_ratio_threshold) & (count <= (max.iter - 1))) {
      SSD_old <- SSD_new
      f0 <- f_new

      T_new <- cbind(rep(1, I_new), t_new)
      # E_new <- matrix(NA, ncol = I_new, nrow = I_new)
      # for (j in 1:I_new) {
      #   E.prepare <- function(t) {
      #     eta_kernel(t - t_new[j, ], gamma)
      #   }
      #   E_new[, j] <- apply(t_new, 1, E.prepare)
      # }

      E_new <- calcE(t_new, gamma)

      # M1_new <- cbind(
      #   2 * E_new %*% W_new %*% E_new + 2 * w_new * E_new,
      #   2 * E_new %*% W_new %*% T_new,
      #   T_new
      # )
      # M2_new <- cbind(
    #   2 * t(T_new) %*% W_new %*% E_new,
      #   2 * t(T_new) %*% W_new %*% T_new,
      #   matrix(0, ncol = d_new + 1, nrow = d_new + 1)
      # )
      # M3_new <- cbind(
      #   t(T_new),
      #   matrix(0, ncol = d_new + 1, nrow = d_new + 1),
      #   matrix(0, ncol = d_new + 1, nrow = d_new + 1)
      # )
      # M_new <- rbind(
      #   M1_new,
      #   M2_new,
      #   M3_new
      # )
      #
      # b_new <- rbind(
      #   2 * E_new %*% W_new %*% X_new,
      #   2 * t(T_new) %*% W_new %*% X_new,
      #   matrix(0, nrow = d_new + 1, ncol = D_new)
      # )
      # sol_new <- ginv(M_new) %*% b_new

      # sol_new <- solve_eq_inter_t(E_new, W_new, T_new, X_new, w_new, d_new, D_new)
      sol_new <- solve_eq(E_new, W_new, T_new, X_new, w_new, d_new, D_new)

      # eta.func <- function(t) {
      #   eta.func.prepare <- function(tau) {
      #     return(eta_kernel(t - tau, gamma))
      #   }
      #   return(matrix(apply(t_new, 1, eta.func.prepare), ncol = 1))
      # }

      f_new <- function(t) {
        return_vec <- as.vector(
          t(sol_new[1:I_new, ]) %*%
            etaFunc(t, t_new, gamma) + t(sol_new[(I_new + 1):(I_new + d_new + 1), ]) %*%
            matrix(c(1, t), ncol = 1)
          )
        return_vec[1] <- t[1]
        return(return_vec)
      }

      t_old <- t_new

      X_initial_guess <- cbind(X_new, t_new)
      projection.index.f0 <- function(x.init) {
        projection_lpme(x.init[1:D_new], f0, x.init[(D_new + 1):(D_new + d_new)])
      }
      t_new <- matrix(
        t(apply(X_initial_guess, 1, projection.index.f0)),
        nrow = I_new
      )
      # t_new <- calc_tnew(X_new, t_new, sol_new, I_new, d_new, gamma)
      # t_new[, 1] <- t_old[, 1]

      X_projection_index <- cbind(X_new, t_new)
      SSD.prepare.again <- function(x.init) {
        return(SSD.prepare(x.init, f_new))
      }

      if (print_plots == TRUE) {
        plot_lpme(df, f_new, t_new, d_new, D_new, time_points)
      }

      # SSD_new <- apply(X_projection_index, 1, SSD.prepare.again) %>%
      #   as.vector() %>%
      #   sum()

      SSD_new <- map(
        1:I_new,
        ~ dist_euclideanC(
          X_new[.x, ],
          fNew(t_new[.x, ], sol_new, t_new, I_new, d_new, gamma)
        ) ^ 2
      ) %>%
        unlist() %>%
        sum()

      SSD_ratio <- abs(SSD_new - SSD_old) / SSD_old
      count <- count + 1

      if (print.SSDs == TRUE) {
        print(
          paste0(
            "SSD = ",
            as.character(round(SSD_new, 4)),
            ", SSD.ratio is ",
            as.character(round(SSD_ratio, 4)),
            ", and this is the ",
            as.character(count),
            "th step of iteration."
          )
        )
      }
    }

    data_initial <- matrix(0, nrow = 1, ncol = D_new + d_new)
    for (i in 1:I_new) {
      index.temp <- which(clusters_full == i)
      length.temp <- length(index.temp)
      X.i <- matrix(df[index.temp, ], nrow = length.temp)
      # X.i <- matrix(df[index.temp, -1], nrow = length.temp)
      t.temp <- matrix(rep(t_new[i, 1], length.temp))
      for (j in 1:d_new) {
        t.temp <- cbind(t.temp, rep(t_new[i, j], length.temp))
      }
      t.temp <- matrix(t.temp[, -1], nrow = length.temp)
      if (length.temp > 0) {
        data_initial <- rbind(data_initial, cbind(X.i, t.temp))
      }
    }
    # for (i in 1:I_new) {
    #   length.temp <- 1
    #   X.i <- matrix(X_new[i, ], nrow = 1)
    #   t_temp <- matrix(rep(t_new[i, 1], length.temp))
    #   for (j in 1:d_new) {
    #     t_temp <- cbind(
    #       t_temp,
    #       rep(t_new[i, j], length.temp)
    #     )
    #   }
    #   t_temp <- matrix(t_temp[, -1], nrow = length.temp)
    #   data_initial <- rbind(data_initial, cbind(X.i, t_temp))
    # }

    data_initial <- data_initial[-1, ]
    proj.para.prepare <- function(data.init) {
      return(
        projection(
          data.init[1:D_new],
          f_new,
          data.init[(D_new + 1):(D_new + d_new)]
        )
      )
    }
    proj_para <- matrix(
      t(apply(data_initial, 1, proj.para.prepare)),
      ncol = d_new
    )
    proj_points <- t(apply(proj_para, 1, f_new))
    diff_data_fit <- apply(
      data_initial[, 1:D_new] - proj_points,
      1,
      norm_euclidean
    )
    MSE_new <- mean(diff_data_fit ^ 2)
    MSE_seq_new[tuning.ind] <- MSE_new
    if (print.MSDs == TRUE) {
      print(
        paste(
          "When gamma = ",
          as.character(w_new),
          ", MSD = ",
          as.character(MSE_new),
          "."
        )
      )
    }

    SOL_new[[tuning.ind]] <- sol_new
    TNEW_new[[tuning.ind]] <- t_new

    if (tuning.ind >= 4) {
      if (
        (MSE_seq_new[tuning.ind] > increase_threshold * MSE_seq_new[tuning.ind - 1]) &
        (MSE_seq_new[tuning.ind - 1] > increase_threshold * MSE_seq_new[tuning.ind - 2]) &
        (MSE_seq_new[tuning.ind - 2] > increase_threshold * MSE_seq_new[tuning.ind - 3])
      ) {
        break
      }
    }
  }

  optimal_ind <- min(which(MSE_seq_new == min(MSE_seq_new)))
  sol_opt <- SOL_new[[optimal_ind]]
  t_new_opt <- TNEW_new[[optimal_ind]]
  # eta.func <- function(t) {
  #   eta.func.prepare <- function(tau) {
  #     return(eta_kernel(t - tau, gamma))
  #   }
  #   return(
  #     matrix(
  #       apply(t_new_opt, 1, eta.func.prepare),
  #       ncol = 1
  #     )
  #   )
  # }
  f.optimal <- function(t) {
    return_vec <- as.vector(
        t(sol_opt[1:I_new, ]) %*%
          etaFunc(t, t_new_opt, gamma) + t(sol_opt[(I_new + 1):(I_new + d_new + 1), ]) %*%
          matrix(c(1, t), ncol = 1)
      )
    return_vec[1] <- t[1]
    return(return_vec)
  }

  if (print.MSDs == TRUE) {
    plot(
      log(tuning.para.seq[1:length(MSE_seq_new)]),
      MSE_seq_new,
      xlab = "Log Gamma",
      ylab = "MSD",
      type = "l"
    )
    lines(
      log(tuning.para.seq[1:length(MSE_seq_new)]),
      MSE_seq_new,
      type = "p",
      pch = 20,
      col = "orange",
      cex = 2
    )
    abline(
      v = log(tuning.para.seq[optimal_ind]),
      lwd = 1.5,
      col = "darkgreen",
      lty = 2
    )
    print(
      paste(
        "The optimal tuning parameter is ",
        as.character(tuning.para.seq[optimal_ind]),
        ", and the MSD of the optimal fit is ",
        as.character(MSE_seq_new[optimal_ind], ".")
      )
    )
  }

  resp <- list(
    embedding_map = f.optimal,
    MSD = MSE_seq_new,
    knots = centers_new,
    weights_of_knots = theta_hat_new,
    clusters = clusters_full,
    coe_kernel = sol_opt[1:I_new, ],
    coe_poly = sol_opt[(I_new + 1):(I_new + d_new + 1), ],
    SOL = SOL_new,
    TNEW = TNEW_new,
    T_parameter = sol_opt,
    Coef = t_new_opt
  )
  return(resp)
}
