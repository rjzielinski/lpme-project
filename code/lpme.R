lpme <- function(df, d, tuning.para.seq = c(0, exp(-15:10)), alpha = 0.05, max.comp = 500, epsilon = 0.05, max.iter = 100, verbose = "all", print_plots = TRUE, SSD_ratio_threshold = 100, increase_threshold = 1.05, init = "full") {
  # df is an N x (D + 1) matrix, with the first column corresponding
  # to the time point at which each observation was collected
  # this matrix should include the observations from all time points

  source("code/pme.R")
  source("code/functions/plot_lpme.R")
  source("code/functions/projection_lpme.R")
  source("code/functions/solve_eq2.R")
  require(plotly)
  require(svMisc)
  require(Rcpp)
  require(nloptr)
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
        init_centers[[idx]] <- init_est_temp$mu
        init_sigma[[idx]] <- init_est_temp$sigma
        init_clusters[[idx]] <- init_est_temp$k.means.result
      }
    }

    init_timevals <- reduce(init_timevals, c)
    init_theta_hat <- reduce(init_theta_hat, c)
    init_centers <- reduce(init_centers, rbind)
    init_sigma <- reduce(init_sigma, c)
    init_W <- diag(init_theta_hat)
    init_X <- init_centers
    init_I <- length(init_theta_hat)

    init_dissimilarity_matrix <- as.matrix(dist(init_X))
    init_isomap <- isomap(init_dissimilarity_matrix, ndim = d, k = 10)
  }

  num_clusters <- rep(0, length(time_points))
  errors <- vector()

  for (idx in 1:length(time_points)) {
    df_temp <- df[df[, 1] == time_points[idx], ]
    if (init %in% c("first", "full")) {
      pme_results[[idx]] <- pme(
        x.obs = df_temp[, -1],
        d = d,
        tuning.para.seq = exp(-20:10),
        initialization = list(
          init_isomap,
          init_theta_hat[init_timevals == time_points[idx]],
          init_centers[init_timevals == time_points[idx], ],
          init_sigma[idx],
          matrix(init_isomap$points[init_timevals == time_points[idx], ], nrow = length(init_theta_hat[init_timevals == time_points[idx]])),
          init_clusters[[idx]]
        ),
        verbose = "none"
      )
    } else {
      pme_results[[idx]] <- pme(
        x.obs = df_temp[, -1],
        d = d,
        verbose = "none"
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

    # r[[idx]] <- cbind(time_points[idx], r_mat)[, -2]
    r[[idx]] <- r_mat
    r2[[idx]] <- time_points[idx]
    errors[idx] <- pme_results[[idx]]$MSD[opt_run]
  }

  r_full <- reduce(r, rbind)

  length_r <- ceiling(max(num_clusters) ^ (1 / d))
  r_bounds <- colMinsMaxs(r_full)
  r_list <- list()
  for (idx in 1:dim(r_bounds)[2]) {
    r_list[[idx]] <- seq(
      r_bounds[1, idx],
      r_bounds[2, idx],
      length.out = length_r
    )
  }
  # r_list[[1]] <- time_points
  r_mat <- as.matrix(expand.grid(r_list))
  n_knots <- nrow(r_mat)
  r_full <- r_mat
  r_full2 <- reduce(r2, rbind)
  r_df <- data.frame(r_full)
  r_df2 <- data.frame(r_full2)
  rnames <- paste0("r", 1:(dim(r_full)[2] - 1))
  # names(r_df) <- c("time", rnames)

  lambda <- vector()
  x_vals <- list()
  spline_coefs <- list()

  gamma <- 4 - d

  for (time_idx in 1:length(time_points)) {
    lambda[time_idx] <- pme_results[[time_idx]]$tuning
    f <- funcs[[time_idx]]
    output <- map(
      1:n_knots,
      ~ f(r_full[.x, ])
    ) %>%
      reduce(rbind)
    x_vals[[time_idx]] <- cbind(time_points[time_idx], output)

    R <- cbind(rep(1, n_knots), r_full)
    E <- calcE(r_full, gamma)

    spline_coefs[[time_idx]] <- solve_eq2(E, R, output, nrow(R), lambda[time_idx], ncol(r_full), ncol(output)) %>%
      t() %>%
      matrix(nrow = 1)
  }

  coef_full <- reduce(spline_coefs, rbind)
  x_test <- reduce(x_vals, rbind)

  # x_test <- reduce(x_test, rbind)
  # x_test_df <- data.frame(x_test)
  # xnames <- paste0("x", 1:(dim(x_test)[2]))
  # names(x_test_df) <- xnames

  # coef_full <- reduce(pme_coefs, rbind)
  # coef_full2 <- reduce(pme_coefs2, rbind)

  # centers_full <- reduce(centers, rbind)
  # centers_full <- cbind(x_test[, 1], centers_full)
  # clusters_full <- unlist(clusters)

  D_coef <- dim(coef_full)[2]
  # D_coef2 <- dim(coef_full2)[2]
  D_new <- dim(x_test)[2] - 1
  D_new2 <- dim(x_test)[2]
  d_new <- dim(r_full)[2]
  d_new2 <- dim(r_full2)[2]
  n <- dim(x_test)[1]
  gamma <- 4 - d_new
  gamma2 <- 4 - d_new2

  # sig_new <- 0.001
  # theta_hat_new <- weight_seq(x_test, x_test, sig_new)
  # theta_hat_new <- init_theta_hat / sum(init_theta_hat)
  # centers_new <- x_test
  # sigma_new <- 0.001
  # W_new <- diag(theta_hat_new)
  # W_new2 <- diag(rep(1 / nrow(r_full2), nrow(r_full2)))
  X_new <- x_test
  # I_new <- length(theta_hat_new)
  I_new <- n

  # dissimilarity_matrix_new <- as.matrix(dist(X_new))
  t_initial <- r_full %>%
    as.matrix()

  MSE_seq_new <- vector()
  SOL_coef <- list()
  TNEW_new <- list()
  coefs <- list()
  x_funs <- list()
  functions <- list()
  func_coef <- list()

  inv_errors <- 1 / errors
  weights <- inv_errors / sum(inv_errors)

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
    W <- diag(weights)
    t_new <- t_initial
    T_new <- cbind(rep(1, n_knots), t_new)
    T_new2 <- cbind(rep(1, nrow(r_full2)), r_full2)

    E_new <- calcE(t_new, gamma)
    E_new2 <- calcE(r_full2, gamma2)

    # sol_coef <- solve_eq2(E_new2, T_new2, coef_full, nrow(coef_full), w_new, d_new2, D_coef)
    sol_coef <- solve_eq(E_new2, W, T_new2, coef_full, w_new, d_new2, D_coef)

    f_coef <- function(t) {
      return_vec <- as.vector(
        # t(sol_coef[1:nrow(coef_full), ]) %*% etaFunc(t, t_new, gamma) +
        t(sol_coef[1:nrow(coef_full), ]) %*% etaFunc(t, r_full2, gamma2) +
          t(sol_coef[(nrow(coef_full) + 1):(nrow(coef_full) + d_new2 + 1), ]) %*% matrix(c(1, t), ncol = 1)
      )
      return(return_vec)
    }

    f_new <- function(t) {
      coefs <- f_coef(t[1])
      coef_mat <- matrix(coefs, n_knots + d_new + 1, byrow = TRUE)
      return_vec <- t(coef_mat[1:n_knots, ]) %*% etaFunc(t[-1], t_new, gamma) +
        t(coef_mat[(n_knots + 1):(n_knots + d_new + 1), ]) %*% matrix(c(1, t[-1]), ncol = 1)
      return(c(t[1], return_vec))
    }

    f0_new <- f_new

    full_t <- expand_grid(time_points, t_new) %>%
      as.matrix(ncol = d + 1)

    X_initial_guess <- cbind(X_new, full_t)
    projection.index.f0 <- function(x.init) {
      projection_lpme(
        x.init[1:D_new2],
        f0_new,
        x.init[(D_new2 + 1):(D_new2 + d_new + 1)]
      )
    }

    t_new2 <- matrix(
      t(apply(X_initial_guess, 1, projection.index.f0)),
      nrow = I_new
    )

    x_fun <- map(
      1:I_new,
      ~ f_new(full_t[.x, ])
    ) %>%
      reduce(rbind)

    # t_new <- calc_tnew(X_new, t_new, sol_new, I_new, d_new, gamma)
    # t_new[, 1] <- t_initial[, 1]

    X_projection_index <- cbind(x_fun, full_t)
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
        f_new(t_new2[.x, ])
      ) ^ 2
    ) %>%
      unlist() %>%
      sum()

    if (print_plots == TRUE) {
      plot_lpme(df, f_new, t_new, d_new, D_new, time_points)
    }

    count <- 1
    SSD_ratio <- 10 * epsilon

    # while ((SSD_ratio > epsilon) & (SSD_ratio <= SSD_ratio_threshold) & (count <= (max.iter - 1))) {
    #
    #   SSD_old <- SSD_new
    #   f0 <- f_new
    #
    #   spline_coefs <- list()
    #
    #   x_vals <- map(
    #     1:I_new,
    #     ~ f_new(full_t[.x, ])
    #   ) %>%
    #     reduce(rbind)
    #
    #   for (time_idx in 1:length(time_points)) {
    #
    #     R <- cbind(rep(1, n_knots), r_mat)
    #     E <- calcE(r_mat, gamma)
    #
    #     spline_coefs[[time_idx]] <- solve_eq2(
    #       E,
    #       R,
    #       x_vals[x_vals[, 1] == time_points[time_idx], -1],
    #       n_knots,
    #       lambda[time_idx],
    #       d_new,
    #       D_new
    #     ) %>%
    #       t() %>%
    #       matrix(nrow = 1)
    #   }
    #
    #   coef_full <- reduce(spline_coefs, rbind)
    #   x_test <- x_vals
    #   X_new <- x_test
    #   # I_new <- length(theta_hat_new)
    #   E_new2 <- calcE(r_full2, gamma2)
    #
    #   sol_coef <- solve_eq2(E_new2, T_new2, coef_full, nrow(coef_full), w_new, d_new2, D_coef)
    #
    #   f_coef <- function(t) {
    #     return_vec <- as.vector(
    #       # t(sol_coef[1:nrow(coef_full), ]) %*% etaFunc(t, t_new, gamma) +
    #       t(sol_coef[1:nrow(coef_full), ]) %*% etaFunc(t, r_full2, gamma2) +
    #         t(sol_coef[(nrow(coef_full) + 1):(nrow(coef_full) + d_new2 + 1), ]) %*% matrix(c(1, t), ncol = 1)
    #     )
    #     return(return_vec)
    #   }
    #
    #   f_new <- function(t) {
    #     coefs <- f_coef(t[1])
    #     coef_mat <- matrix(coefs, n_knots + d_new + 1, byrow = TRUE)
    #     return_vec <- t(coef_mat[1:n_knots, ]) %*% etaFunc(t[-1], t_new, gamma) +
    #       t(coef_mat[(n_knots + 1):(n_knots + d_new + 1), ]) %*% matrix(c(1, t[-1]), ncol = 1)
    #     return(c(t[1], return_vec))
    #   }
    #
    #   f0_new <- f_new
    #
    #   X_initial_guess <- cbind(X_new, full_t)
    #   projection.index.f0 <- function(x.init) {
    #     projection_lpme(
    #       x.init[1:D_new2],
    #       f0_new,
    #       x.init[(D_new2 + 1):(D_new2 + d_new + 1)]
    #     )
    #   }
    #
    #   t_old <- t_new2
    #
    #   t_new2 <- matrix(
    #     t(apply(X_initial_guess, 1, projection.index.f0)),
    #     nrow = I_new
    #   )
    #
    #   SSD_new <- map(
    #     1:I_new,
    #     ~ dist_euclideanC(
    #       X_new[.x, ],
    #       # fNew(t_new[.x, ], sol_new, t_new, I_new, d_new, gamma)
    #       f_new(t_new2[.x, ])
    #     ) ^ 2
    #   ) %>%
    #     unlist() %>%
    #     sum()
    #
    #   if (print_plots == TRUE) {
    #     plot_lpme(df, f_new, r_mat, d_new, D_new, time_points)
    #   }
    #
    #   SSD_ratio <- abs(SSD_new - SSD_old) / SSD_old
    #   count <- count + 1
    #
    #   if (print.SSDs == TRUE) {
    #     print(
    #       paste0(
    #         "SSD = ",
    #         as.character(round(SSD_new, 4)),
    #         ", SSD.ratio is ",
    #         as.character(round(SSD_ratio, 4)),
    #         ", and this is the ",
    #         as.character(count),
    #         "th step of iteration."
    #       )
    #     )
    #   }
    # }

    nearest_x <- calc_nearest_x(df, x_test)
    init_param <- calc_init_param(df, t_new2, nearest_x)

    data_initial <- cbind(df, init_param)

    cv_mse <- vector()
    for (time_idx in 1:length(time_points)) {
      r_bounds <- colMinsMaxs(t_new2[, -1])
      r_list <- list()
      for (idx in 1:dim(r_bounds)[2]) {
        r_list[[idx]] <- seq(
          r_bounds[1, idx],
          r_bounds[2, idx],
          length.out = n_knots ^ (1 / d)
        )
      }
      # r_list[[1]] <- time_points
      r_mat <- as.matrix(expand.grid(r_list))
      r_full_cv <- expand_grid(time_points[-time_idx], r_mat) %>%
        as.matrix()

      spline_coefs <- list()

      x_vals <- map(
        1:nrow(r_full_cv),
        ~ f_new(r_full_cv[.x, ])
      ) %>%
        reduce(rbind)

      for (idx in 1:length(time_points)) {
        if (time_idx == idx) {
          next
        }

        R <- cbind(rep(1, n_knots), r_mat)
        E <- calcE(r_mat, gamma)

        spline_coefs[[idx]] <- solve_eq2(
          E,
          R,
          x_vals[x_vals[, 1] == time_points[idx], -1],
          n_knots,
          lambda[idx],
          d_new,
          D_new
        ) %>%
          t() %>%
          matrix(nrow = 1)
      }

      coef_cv <- reduce(spline_coefs, rbind)
      x_test_cv <- x_vals
      X_new_cv <- x_test_cv
      # I_new <- length(theta_hat_new)
      I_new_cv <- nrow(X_new_cv)

      r_cv <- matrix(
        r_full2[r_full2[, 1] != time_points[time_idx], ],
        ncol = 1
      )

      E_new_cv <- calcE(r_cv, gamma2)
      T_new_cv <- T_new2[T_new2[, 2] != time_points[time_idx], ]

      sol_coef_cv <- solve_eq2(E_new_cv, T_new_cv, coef_cv, nrow(coef_cv), w_new, d_new2, D_coef)

      f_coef_cv <- function(t) {
        return_vec <- as.vector(
          # t(sol_coef[1:nrow(coef_full), ]) %*% etaFunc(t, t_new, gamma) +
          t(sol_coef_cv[1:nrow(coef_cv), ]) %*% etaFunc(t, r_cv, gamma2) +
            t(sol_coef[(nrow(coef_cv) + 1):(nrow(coef_cv) + d_new2 + 1), ]) %*% matrix(c(1, t), ncol = 1)
        )
        return(return_vec)
      }

      f_new_cv <- function(t) {
        coefs <- f_coef_cv(t[1])
        coef_mat <- matrix(coefs, n_knots + d_new + 1, byrow = TRUE)
        return_vec <- t(coef_mat[1:n_knots, ]) %*% etaFunc(t[-1], t_new, gamma) +
          t(coef_mat[(n_knots + 1):(n_knots + d_new + 1), ]) %*% matrix(c(1, t[-1]), ncol = 1)
        return(c(t[1], return_vec))
      }

      temp_data_initial <- data_initial[data_initial[, 1] == time_points[time_idx], ]
      proj_para_cv <- map(
        1:nrow(temp_data_initial),
        ~ try({
          projection_lpme(
            temp_data_initial[.x, 1:D_new2],
            f_new_cv,
            temp_data_initial[.x, (D_new2 + 1):(D_new2 + d_new + 1)]
          ) %>%
            t()
        })
      )

      errors <- vector()
      for (i in 1:length(proj_para_cv)) {
        errors[i] <- sum(class(proj_para_cv[[i]]) == "try-error") > 0
      }

      if (sum(errors) > 0) {
        proj_para_cv <- reduce(proj_para_cv, rbind)[-which(errors), ]
      } else {
        proj_para_cv <- reduce(proj_para_cv, rbind)
      }

      proj_points_cv <- t(apply(proj_para_cv, 1, f_new_cv))

      if (sum(errors) > 0) {
        temp_data_initial <- temp_data_initial[-which(errors), ]
      }

      proj_error_cv <- dist_euclideanC_vec(
        temp_data_initial[, 1:D_new2],
        proj_points_cv
      ) %>%
        as.vector()
      cv_mse[time_idx] <- mean(proj_error_cv ^ 2)
    }

    df_n <- sapply(time_points, function(x) nrow(df[df[, 1] == x, ]))
    MSE_new <- weighted.mean(cv_mse, df_n)
    MSE_seq_new[tuning.ind] <- MSE_new


    # proj_para <- map(
    #   1:nrow(data_initial),
    #   ~ projection_lpme(
    #     data_initial[.x, 1:D_new2],
    #     f_new,
    #     data_initial[.x, (D_new2 + 1):(D_new2 + d_new + 1)]
    #   ) %>%
    #     t()
    # ) %>%
    #   reduce(rbind)
    #
    # proj_points <- t(apply(proj_para, 1, f_new))
    #
    # proj_error <- map(
    #   1:nrow(df),
    #   ~ dist_euclideanC(df[.x, ], proj_points[.x, ])
    # ) %>%
    #   reduce(c)
    #
    # MSE_new <- mean(proj_error ^ 2)
    # MSE_seq_new[tuning.ind] <- MSE_new
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

    SOL_coef[[tuning.ind]] <- sol_coef
    TNEW_new[[tuning.ind]] <- t_new2
    coefs[[tuning.ind]] <- r_mat
    x_funs[[tuning.ind]] <- x_fun
    functions[[tuning.ind]] <- f_new
    func_coef[[tuning.ind]] <- f_coef

    # if (tuning.ind >= 4) {
    #   if (
    #     (MSE_seq_new[tuning.ind] > increase_threshold * MSE_seq_new[tuning.ind - 1]) &
    #     (MSE_seq_new[tuning.ind - 1] > increase_threshold * MSE_seq_new[tuning.ind - 2]) &
    #     (MSE_seq_new[tuning.ind - 2] > increase_threshold * MSE_seq_new[tuning.ind - 3])
    #   ) {
    #     break
    #   }
    # }
  }

  optimal_ind <- min(which(MSE_seq_new == min(MSE_seq_new)))
  sol_opt <- SOL_coef[[optimal_ind]]
  t_new_opt <- TNEW_new[[optimal_ind]]
  coefs_opt <- coefs[[optimal_ind]]
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
  # f.optimal <- function(t) {
  #   return_vec <- as.vector(
  #       t(sol_opt[1:I_new, ]) %*%
  #         etaFunc(t, t_new_opt, gamma) + t(sol_opt[(I_new + 1):(I_new + d_new + 1), ]) %*%
  #         matrix(c(1, t), ncol = 1)
  #     )
  #   return_vec[1] <- t[1]
  #   return(return_vec)
  # }

  f.optimal <- functions[[optimal_ind]]

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
    # knots = centers_new,
    # weights_of_knots = theta_hat_new,
    # clusters = clusters_full,
    # coe_kernel = sol_opt,
    # coe_poly = sol_opt2,
    SOL = SOL_coef,
    TNEW = TNEW_new,
    coefs = coefs,
    x_vals = x_funs,
    T_parameter = sol_opt,
    Coef = t_new_opt,
    f_out = functions,
    f_coef = func_coef
  )
  return(resp)
}

