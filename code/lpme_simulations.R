library(tidyverse)
library(plotly)
library(pracma)
library(foreach)
library(doParallel)
library(furrr)
library(progress)
library(pme)
library(princurve)

source("code/functions/sim_data.R")
source("code/functions/calc_pme_est.R")
source("code/functions/calc_lpme_est.R")

### Simulation Case 1

sim_error_case1 <- function(max_time, interval, amp_noise, shape_noise, time_change, time_trend, n, run = 1, print_plots = FALSE) {
  require(tidyverse)
  require(pme)
  require(princurve)
  source("code/functions/sim_data.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  time_vals <- seq(0, max_time, interval)
  sim_list <- lapply(
    time_vals,
    sim_data,
    case = 1,
    noise = 0.15,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    N = n,
    time_change = time_change,
    time_trend = time_trend
  )

  sim_df <- matrix(ncol = ncol(sim_list[[1]][[1]]))
  true_vals <- matrix(ncol = ncol(sim_list[[1]][[2]]))
  for (i in 1:length(sim_list)) {
    sim_df <- rbind(sim_df, sim_list[[i]][[1]])
    true_vals <- rbind(true_vals, sim_list[[i]][[2]])
  }
  sim_df <- sim_df[-1, ]
  true_vals <- true_vals[-1, ]

  lpme_result <- lpme(
    sim_df,
    1,
    smoothing_method = "spline",
    print_plots = print_plots,
    verbose = TRUE
  )
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)

  # lpme_result_gp <- lpme(
  #   sim_df,
  #   1,
  #   smoothing_method = "gp",
  #   print_plots = print_plots,
  #   verbose = TRUE
  # )
  # lpme_vals_gp <- calc_lpme_est(lpme_result_gp, sim_df)

  pme_result <- list()
  pme_vals <- list()
  smoothing_options <- c("smooth_spline", "lowess", "periodic_lowess")
  principal_curve_result <- list()
  principal_curve_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = FALSE)
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
    curves <- list()
    curve_vals <- list()
    curve_error <- vector()
    for (smoother_idx in 1:length(smoothing_options)) {
      curves[[smoother_idx]] <- principal_curve(
        temp_data,
        smoother = smoothing_options[smoother_idx]
      )
      curve_vals[[smoother_idx]] <- cbind(time_vals[t], curves[[smoother_idx]]$s)
      curve_error[smoother_idx] <- map(
        1:nrow(true_vals[true_vals[, 1] == time_vals[t], ]),
        ~ dist_euclidean(true_vals[true_vals[, 1] == time_vals[t], ][.x, ], curve_vals[[smoother_idx]][.x, ])^2
      ) %>%
        unlist() %>%
        mean()
    }
    opt_curve <- which.min(curve_error)
    principal_curve_result[[t]] <- curves[[opt_curve]]
    principal_curve_vals[[t]] <- curve_vals[[opt_curve]]
  }
  pme_vals <- reduce(pme_vals, rbind)
  principal_curve_vals <- reduce(principal_curve_vals, rbind)

  p <- plot_ly(
    x = lpme_vals[, 2],
    y = lpme_vals[, 3],
    z = lpme_vals[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  ) %>%
    add_markers(
      x = pme_vals[, 2],
      y = pme_vals[, 3],
      z = pme_vals[, 1]
    ) %>%
    add_markers(
      x = principal_curve_vals[, 2],
      y = principal_curve_vals[, 3],
      z = principal_curve_vals[, 1]
    ) %>%
    add_markers(
      x = true_vals[, 2],
      y = true_vals[, 3],
      z = true_vals[, 1]
    )

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], pme_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  principal_curve_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], principal_curve_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], lpme_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  # lpme_error_gp <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, ], lpme_vals_gp[.x, ])^2
  # ) %>%
  #   unlist() %>%
  #   mean()

  sim_case1 <- list(
    df = sim_df,
    times = time_vals,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    n = n,
    lpme_result = lpme_result,
    # lpme_gp_result = lpme_result_gp,
    pme_results = pme_result,
    principal_curve_results = principal_curve_result,
    lpme_error = lpme_error,
    # lpme_gp_error = lpme_error_gp,
    pme_error = pme_error,
    principal_curve_error = principal_curve_error,
    plot = p
  )
  sim_dir <- "simulations/case1/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, recursive = TRUE)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_ampnoise_",
    str_pad(as.character(100 * amp_noise), 3, side = "left", pad = "0"),
    "_pernoise_",
    str_pad(as.character(100 * shape_noise), 3, side = "left", pad = "0"),
    "_n_",
    str_pad(as.character(n), 4, side = "left", pad = "0"),
    "_",
    time_trend,
    str_pad(as.character(100 * time_change), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case1,
    paste0(sim_dir, filename)
  )
  return(0)
}

### Simulation Case 2

sim_error_case2 <- function(max_time, interval, amp_noise, shape_noise, time_change, time_trend, n, run = 1, print_plots = FALSE) {
  require(tidyverse)
  require(pme)
  require(princurve)
  source("code/functions/sim_data.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")

  time_vals <- seq(0, max_time, interval)
  sim_list <- lapply(
    time_vals,
    sim_data,
    case = 2,
    noise = 0.15,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    N = n,
    time_change = time_change,
    time_trend = time_trend
  )

  sim_df <- matrix(ncol = ncol(sim_list[[1]][[1]]))
  true_vals <- matrix(ncol = ncol(sim_list[[1]][[2]]))
  for (i in 1:length(sim_list)) {
    sim_df <- rbind(sim_df, sim_list[[i]][[1]])
    true_vals <- rbind(true_vals, sim_list[[i]][[2]])
  }
  sim_df <- sim_df[-1, ]
  true_vals <- true_vals[-1, ]


  lpme_result <- lpme(
    sim_df,
    1,
    smoothing_method = "spline",
    print_plots = print_plots,
    verbose = TRUE
  )
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  # lpme_result_gp <- lpme(
  #   sim_df,
  #   1,
  #   smoothing_method = "gp",
  #   print_plots = print_plots,
  #   verbose = TRUE
  # )
  # lpme_vals_gp <- calc_lpme_est(lpme_result_gp, sim_df)
  pme_result <- list()
  pme_vals <- list()
  smoothing_options <- c("smooth_spline", "lowess", "periodic_lowess")
  principal_curve_result <- list()
  principal_curve_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
    curves <- list()
    curve_vals <- list()
    curve_error <- vector()
    for (smoother_idx in 1:length(smoothing_options)) {
      curves[[smoother_idx]] <- principal_curve(
        temp_data,
        smoother = smoothing_options[smoother_idx]
      )
      curve_vals[[smoother_idx]] <- cbind(time_vals[t], curves[[smoother_idx]]$s)
      curve_error[smoother_idx] <- map(
        1:nrow(true_vals[true_vals[, 1] == time_vals[t], ]),
        ~ dist_euclidean(true_vals[true_vals[, 1] == time_vals[t], ][.x, ], curve_vals[[smoother_idx]][.x, ])^2
      ) %>%
        unlist() %>%
        mean()
    }
    opt_curve <- which.min(curve_error)
    principal_curve_result[[t]] <- curves[[opt_curve]]
    principal_curve_vals[[t]] <- curve_vals[[opt_curve]]
  }
  pme_vals <- reduce(pme_vals, rbind)
  principal_curve_vals <- reduce(principal_curve_vals, rbind)

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], pme_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], lpme_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  # lpme_error_gp <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, ], lpme_vals_gp[.x, ])^2
  # ) %>%
  #   unlist() %>%
  #   mean()

  principal_curve_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], principal_curve_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  p <- plot_ly(
    x = lpme_vals[, 2],
    y = lpme_vals[, 3],
    z = lpme_vals[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  ) %>%
    add_markers(
      x = pme_vals[, 2],
      y = pme_vals[, 3],
      z = pme_vals[, 1]
    ) %>%
    add_markers(
      x = principal_curve_vals[, 2],
      y = principal_curve_vals[, 3],
      z = principal_curve_vals[, 1]
    ) %>%
    add_markers(
      x = true_vals[, 2],
      y = true_vals[, 3],
      z = true_vals[, 1]
    )

  sim_case2 <- list(
    df = sim_df,
    times = time_vals,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    n = n,
    lpme_result = lpme_result,
    # lpme_gp_result = lpme_result_gp,
    pme_results = pme_result,
    principal_curve_results = principal_curve_result,
    lpme_error = lpme_error,
    # lpme_gp_error = lpme_error_gp,
    pme_error = pme_error,
    principal_curve_error = principal_curve_error,
    plot = p
  )
  sim_dir <- "simulations/case2/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, recursive = TRUE)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_ampnoise_",
    str_pad(as.character(100 * amp_noise), 3, side = "left", pad = "0"),
    "_pernoise_",
    str_pad(as.character(100 * shape_noise), 3, side = "left", pad = "0"),
    "_n_",
    str_pad(as.character(n), 4, side = "left", pad = "0"),
    "_",
    time_trend,
    str_pad(as.character(100 * time_change), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case2,
    paste0(sim_dir, filename)
  )
  return(0)
}

### Simulation Case 3

sim_error_case3 <- function(max_time, interval, amp_noise, shape_noise, time_change, time_trend, n, run = 1, print_plots = FALSE) {
  require(tidyverse)
  require(pme)
  require(princurve)
  source("code/functions/sim_data.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")

  time_vals <- seq(0, max_time, interval)
  sim_list <- lapply(
    time_vals,
    sim_data,
    case = 3,
    noise = 0.15,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    N = n,
    time_change = time_change,
    time_trend = time_trend
  )

  sim_df <- matrix(ncol = ncol(sim_list[[1]][[1]]))
  true_vals <- matrix(ncol = ncol(sim_list[[1]][[2]]))
  for (i in 1:length(sim_list)) {
    sim_df <- rbind(sim_df, sim_list[[i]][[1]])
    true_vals <- rbind(true_vals, sim_list[[i]][[2]])
  }
  sim_df <- sim_df[-1, ]
  true_vals <- true_vals[-1, ]


  pol <- sim_df[, -1] %>%
    cart2pol()
  sim_df <- cbind(sim_df, pol)[, -5]

  lpme_result <- lpme(
    sim_df,
    1,
    smoothing_method = "spline",
    print_plots = print_plots,
    verbose = TRUE
  )
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  # lpme_result_gp <- lpme(
  #   sim_df,
  #   1,
  #   smoothing_method = "gp",
  #   print_plots = print_plots,
  #   verbose = TRUE
  # )
  # lpme_vals_gp <- calc_lpme_est(lpme_result_gp, sim_df)
  pme_result <- list()
  pme_vals <- list()
  smoothing_options <- c("smooth_spline", "lowess", "periodic_lowess")
  principal_curve_result <- list()
  principal_curve_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
    curves <- list()
    curve_vals <- list()
    curve_error <- vector()
    for (smoother_idx in 1:length(smoothing_options)) {
      curves[[smoother_idx]] <- principal_curve(
        temp_data,
        smoother = smoothing_options[smoother_idx]
      )
      curve_vals[[smoother_idx]] <- cbind(time_vals[t], curves[[smoother_idx]]$s)
      curve_error[smoother_idx] <- map(
        1:nrow(true_vals[true_vals[, 1] == time_vals[t], ]),
        ~ dist_euclidean(true_vals[true_vals[, 1] == time_vals[t], ][.x, ], curve_vals[[smoother_idx]][.x, ])^2
      ) %>%
        unlist() %>%
        mean()
    }
    opt_curve <- which.min(curve_error)
    principal_curve_result[[t]] <- curves[[opt_curve]]
    principal_curve_vals[[t]] <- curve_vals[[opt_curve]]
  }
  pme_vals <- reduce(pme_vals, rbind)
  principal_curve_vals <- reduce(principal_curve_vals, rbind)

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], pme_vals[.x, 1:3])^2
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], lpme_vals[.x, 1:3])^2
  ) %>%
    unlist() %>%
    mean()

  # lpme_error_gp <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, ], lpme_vals_gp[.x, 1:3])^2
  # ) %>%
  #   unlist() %>%
  #   mean()

  principal_curve_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], principal_curve_vals[.x, 1:3])^2
  ) %>%
    unlist() %>%
    mean()

  p <- plot_ly(
    x = lpme_vals[, 2],
    y = lpme_vals[, 3],
    z = lpme_vals[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  ) %>%
    add_markers(
      x = pme_vals[, 2],
      y = pme_vals[, 3],
      z = pme_vals[, 1]
    ) %>%
    add_markers(
      x = principal_curve_vals[, 2],
      y = principal_curve_vals[, 3],
      z = principal_curve_vals[, 1]
    ) %>%
    add_markers(
      x = true_vals[, 2],
      y = true_vals[, 3],
      z = true_vals[, 1]
    )


  sim_case3 <- list(
    df = sim_df,
    times = time_vals,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    n = n,
    lpme_result = lpme_result,
    # lpme_gp_result = lpme_result_gp,
    pme_results = pme_result,
    principal_curve_results = principal_curve_result,
    lpme_error = lpme_error,
    # lpme_gp_error = lpme_error_gp,
    pme_error = pme_error,
    principal_curve_error = principal_curve_error,
    plot = p
  )
  sim_dir <- "simulations/case3/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, recursive = TRUE)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_ampnoise_",
    str_pad(as.character(100 * amp_noise), 3, side = "left", pad = "0"),
    "_pernoise_",
    str_pad(as.character(100 * shape_noise), 3, side = "left", pad = "0"),
    "_n_",
    str_pad(as.character(n), 4, side = "left", pad = "0"),
    "_",
    time_trend,
    str_pad(as.character(100 * time_change), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case3,
    paste0(sim_dir, filename)
  )
  return(0)
}

### Simulation Case 4

sim_error_case4 <- function(max_time, interval, amp_noise, shape_noise, time_change, time_trend, n, run = 1, print_plots = FALSE) {
  require(tidyverse)
  require(pme)
  require(princurve)
  source("code/functions/sim_data.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")

  time_vals <- seq(0, max_time, interval)
  sim_list <- lapply(
    time_vals,
    sim_data,
    case = 4,
    noise = 0.15,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    N = n,
    time_change = time_change,
    time_trend = time_trend
  )

  sim_df <- matrix(ncol = ncol(sim_list[[1]][[1]]))
  true_vals <- matrix(ncol = ncol(sim_list[[1]][[2]]))
  for (i in 1:length(sim_list)) {
    sim_df <- rbind(sim_df, sim_list[[i]][[1]])
    true_vals <- rbind(true_vals, sim_list[[i]][[2]])
  }
  sim_df <- sim_df[-1, ]
  true_vals <- true_vals[-1, ]


  pol <- sim_df[, -1] %>%
    cart2pol()
  # pol[, 1] <- ifelse(pol[, 1] < 0, pol[, 1] + (2 * pi), pol[, 1])
  min_theta <- map(
    time_vals,
    ~ min(pol[(sim_df[, 1] == .x) & pol[, 1] > 0, 1])
  ) %>%
    reduce(c)
  rotation <- map(
    1:nrow(sim_df),
    ~ - ((5/6) * pi) - min_theta[which(time_vals == sim_df[.x, 1])]
  ) %>%
    reduce(c)
  pol[, 1] <- pol[, 1] + rotation
  sim_df2 <- pol2cart(pol)
  sim_df <- cbind(sim_df, cart2pol(sim_df2))

  lpme_result <- lpme(
    sim_df,
    1,
    smoothing_method = "spline",
    print_plots = print_plots,
    verbose = TRUE
  )
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  # lpme_result_gp <- lpme(
  #   sim_df,
  #   1,
  #   smoothing_method = "gp",
  #   print_plots = print_plots,
  #   verbose = TRUE
  # )
  # lpme_vals_gp <- calc_lpme_est(lpme_result_gp, sim_df)
  pme_result <- list()
  pme_vals <- list()
  smoothing_options <- c("smooth_spline", "lowess", "periodic_lowess")
  principal_curve_result <- list()
  principal_curve_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
    curves <- list()
    curve_vals <- list()
    curve_error <- vector()
    for (smoother_idx in 1:length(smoothing_options)) {
      curves[[smoother_idx]] <- principal_curve(
        temp_data,
        smoother = smoothing_options[smoother_idx]
      )
      curve_vals[[smoother_idx]] <- cbind(time_vals[t], curves[[smoother_idx]]$s)
      curve_error[smoother_idx] <- map(
        1:nrow(true_vals[true_vals[, 1] == time_vals[t], ]),
        ~ dist_euclidean(true_vals[true_vals[, 1] == time_vals[t], ][.x, ], curve_vals[[smoother_idx]][.x, ])^2
      ) %>%
        unlist() %>%
        mean()
    }
    opt_curve <- which.min(curve_error)
    principal_curve_result[[t]] <- curves[[opt_curve]]
    principal_curve_vals[[t]] <- curve_vals[[opt_curve]]
  }
  pme_vals <- reduce(pme_vals, rbind)
  principal_curve_vals <- reduce(principal_curve_vals, rbind)

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], pme_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], lpme_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  # lpme_error_gp <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, ], lpme_vals_gp[.x, ])^2
  # ) %>%
  #   unlist() %>%
  #   mean()

  principal_curve_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], principal_curve_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  p <- plot_ly(
    x = lpme_vals[, 2],
    y = lpme_vals[, 3],
    z = lpme_vals[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  ) %>%
    add_markers(
      x = pme_vals[, 2],
      y = pme_vals[, 3],
      z = pme_vals[, 1]
    ) %>%
    add_markers(
      x = principal_curve_vals[, 2],
      y = principal_curve_vals[, 3],
      z = principal_curve_vals[, 1]
    ) %>%
    add_markers(
      x = true_vals[, 2],
      y = true_vals[, 3],
      z = true_vals[, 1]
    )


  sim_case4 <- list(
    df = sim_df,
    times = time_vals,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    n = n,
    lpme_result = lpme_result,
    # lpme_gp_result = lpme_result_gp,
    pme_results = pme_result,
    principal_curve_results = principal_curve_result,
    lpme_error = lpme_error,
    # lpme_gp_error = lpme_error_gp,
    pme_error = pme_error,
    principal_curve_error = principal_curve_error,
    plot = p
  )
  sim_dir <- "simulations/case4/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, recursive = TRUE)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_ampnoise_",
    str_pad(as.character(100 * amp_noise), 3, side = "left", pad = "0"),
    "_pernoise_",
    str_pad(as.character(100 * shape_noise), 3, side = "left", pad = "0"),
    "_n_",
    str_pad(as.character(n), 4, side = "left", pad = "0"),
    "_",
    time_trend,
    str_pad(as.character(100 * time_change), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case4,
    paste0(sim_dir, filename)
  )
  return(0)
}

### Simulation Case 5

sim_error_case5 <- function(max_time, interval, amp_noise, shape_noise, time_change, time_trend, n, run = 1, print_plots = FALSE) {
  require(tidyverse)
  require(pme)
  require(princurve)
  source("code/functions/sim_data.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")

  time_vals <- seq(0, max_time, interval)
  sim_list <- lapply(
    time_vals,
    sim_data,
    case = 5,
    noise = 0.15,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    N = n,
    time_change = time_change,
    time_trend = time_trend
  )

  sim_df <- matrix(ncol = ncol(sim_list[[1]][[1]]))
  true_vals <- matrix(ncol = ncol(sim_list[[1]][[2]]))
  for (i in 1:length(sim_list)) {
    sim_df <- rbind(sim_df, sim_list[[i]][[1]])
    true_vals <- rbind(true_vals, sim_list[[i]][[2]])
  }
  sim_df <- sim_df[-1, ]
  true_vals <- true_vals[-1, ]


  lpme_result <- lpme(
    sim_df,
    1,
    smoothing_method = "spline",
    print_plots = print_plots,
    verbose = TRUE
  )
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)

  # lpme_result_gp <- lpme(
  #   sim_df,
  #   1,
  #   smoothing_method = "gp",
  #   print_plots = print_plots,
  #   verbose = TRUE
  # )
  # lpme_vals_gp <- calc_lpme_est(lpme_result_gp, sim_df)
  pme_result <- list()
  pme_vals <- list()
  smoothing_options <- c("smooth_spline", "lowess", "periodic_lowess")
  principal_curve_result <- list()
  principal_curve_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
    curves <- list()
    curve_vals <- list()
    curve_error <- vector()
    for (smoother_idx in 1:length(smoothing_options)) {
      curves[[smoother_idx]] <- principal_curve(
        temp_data,
        smoother = smoothing_options[smoother_idx]
      )
      curve_vals[[smoother_idx]] <- cbind(time_vals[t], curves[[smoother_idx]]$s)
      curve_error[smoother_idx] <- map(
        1:nrow(true_vals[true_vals[, 1] == time_vals[t], ]),
        ~ dist_euclidean(true_vals[true_vals[, 1] == time_vals[t], ][.x, ], curve_vals[[smoother_idx]][.x, ])^2
      ) %>%
        unlist() %>%
        mean()
    }
    opt_curve <- which.min(curve_error)
    principal_curve_result[[t]] <- curves[[opt_curve]]
    principal_curve_vals[[t]] <- curve_vals[[opt_curve]]
  }
  pme_vals <- reduce(pme_vals, rbind)
  principal_curve_vals <- reduce(principal_curve_vals, rbind)

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], pme_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], lpme_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  # lpme_error_gp <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, ], lpme_vals_gp[.x, ])^2
  # ) %>%
  #   unlist() %>%
  #   mean()

  principal_curve_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], principal_curve_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  p <- plot_ly(
    x = lpme_vals[, 2],
    y = lpme_vals[, 3],
    z = lpme_vals[, 4],
    frame = lpme_vals[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  ) %>%
    add_markers(
      x = pme_vals[, 2],
      y = pme_vals[, 3],
      z = pme_vals[, 4],
      frame = pme_vals[, 1]
    ) %>%
    add_markers(
      x = principal_curve_vals[, 2],
      y = principal_curve_vals[, 3],
      z = principal_curve_vals[, 4],
      frame = principal_curve_vals[, 1]
    ) %>%
    add_markers(
      x = true_vals[, 2],
      y = true_vals[, 3],
      z = true_vals[, 4],
      frame = true_vals[, 1]
    )

  sim_case5 <- list(
    df = sim_df,
    times = time_vals,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    n = n,
    lpme_result = lpme_result,
    # lpme_gp_result = lpme_result_gp,
    pme_results = pme_result,
    principal_curve_results = principal_curve_result,
    lpme_error = lpme_error,
    # lpme_gp_error = lpme_error_gp,
    pme_error = pme_error,
    principal_curve_error = principal_curve_error,
    plot = p
  )
  sim_dir <- "simulations/case5/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, recursive = TRUE)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_ampnoise_",
    str_pad(as.character(100 * amp_noise), 3, side = "left", pad = "0"),
    "_pernoise_",
    str_pad(as.character(100 * shape_noise), 3, side = "left", pad = "0"),
    "_n_",
    str_pad(as.character(n), 4, side = "left", pad = "0"),
    "_",
    time_trend,
    str_pad(as.character(100 * time_change), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case5,
    paste0(sim_dir, filename)
  )
  return(0)
}

### Simulation Case 6

sim_error_case6 <- function(max_time, interval, amp_noise, shape_noise, time_change, time_trend, n, run = 1, print_plots = FALSE) {
  require(tidyverse)
  require(pme)
  require(princurve)
  source("code/functions/sim_data.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")

  time_vals <- seq(0, max_time, interval)
  sim_list <- lapply(
    time_vals,
    sim_data,
    case = 6,
    noise = 0.15,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    N = n,
    time_change = time_change,
    time_trend = time_trend
  )

  sim_df <- matrix(ncol = ncol(sim_list[[1]][[1]]))
  true_vals <- matrix(ncol = ncol(sim_list[[1]][[2]]))
  for (i in 1:length(sim_list)) {
    sim_df <- rbind(sim_df, sim_list[[i]][[1]])
    true_vals <- rbind(true_vals, sim_list[[i]][[2]])
  }
  sim_df <- sim_df[-1, ]
  true_vals <- true_vals[-1, ]


  lpme_result <- lpme(
    sim_df,
    1,
    smoothing_method = "spline",
    print_plots = print_plots,
    verbose = TRUE
  )
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)

  # lpme_result_gp <- lpme(
  #   sim_df,
  #   1,
  #   smoothing_method = "gp",
  #   print_plots = print_plots,
  #   verbose = TRUE
  # )
  # lpme_vals_gp <- calc_lpme_est(lpme_result_gp, sim_df)

  pme_result <- list()
  pme_vals <- list()
  smoothing_options <- c("smooth_spline", "lowess", "periodic_lowess")
  principal_curve_result <- list()
  principal_curve_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
    curves <- list()
    curve_vals <- list()
    curve_error <- vector()
    for (smoother_idx in 1:length(smoothing_options)) {
      curves[[smoother_idx]] <- principal_curve(
        temp_data,
        smoother = smoothing_options[smoother_idx]
      )
      curve_vals[[smoother_idx]] <- cbind(time_vals[t], curves[[smoother_idx]]$s)
      curve_error[smoother_idx] <- map(
        1:nrow(true_vals[true_vals[, 1] == time_vals[t], ]),
        ~ dist_euclidean(true_vals[true_vals[, 1] == time_vals[t], ][.x, ], curve_vals[[smoother_idx]][.x, ])^2
      ) %>%
        unlist() %>%
        mean()
    }
    opt_curve <- which.min(curve_error)
    principal_curve_result[[t]] <- curves[[opt_curve]]
    principal_curve_vals[[t]] <- curve_vals[[opt_curve]]
  }
  pme_vals <- reduce(pme_vals, rbind)
  principal_curve_vals <- reduce(principal_curve_vals, rbind)

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], pme_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], lpme_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  # lpme_error_gp <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, ], lpme_vals_gp[.x, ])^2
  # ) %>%
  #   unlist() %>%
  #   mean()

  principal_curve_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], principal_curve_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  p <- plot_ly(
    x = lpme_vals[, 2],
    y = lpme_vals[, 3],
    z = lpme_vals[, 4],
    frame = lpme_vals[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  ) %>%
    add_markers(
      x = pme_vals[, 2],
      y = pme_vals[, 3],
      z = pme_vals[, 4],
      frame = pme_vals[, 1]
    ) %>%
    add_markers(
      x = principal_curve_vals[, 2],
      y = principal_curve_vals[, 3],
      z = principal_curve_vals[, 4],
      frame = principal_curve_vals[, 1]
    ) %>%
    add_markers(
      x = true_vals[, 2],
      y = true_vals[, 3],
      z = true_vals[, 4],
      frame = true_vals[, 1]
    )

  sim_case6 <- list(
    df = sim_df,
    times = time_vals,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    n = n,
    lpme_result = lpme_result,
    # lpme_gp_result = lpme_result_gp,
    pme_results = pme_result,
    principal_curve_results = principal_curve_result,
    lpme_error = lpme_error,
    # lpme_gp_error = lpme_error_gp,
    pme_error = pme_error,
    principal_curve_error = principal_curve_error,
    plot = p
  )

  sim_dir <- "simulations/case6/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, recursive = TRUE)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_ampnoise_",
    str_pad(as.character(100 * amp_noise), 3, side = "left", pad = "0"),
    "_pernoise_",
    str_pad(as.character(100 * shape_noise), 3, side = "left", pad = "0"),
    "_n_",
    str_pad(as.character(n), 4, side = "left", pad = "0"),
    "_",
    time_trend,
    str_pad(as.character(100 * time_change), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )

  saveRDS(
    sim_case6,
    paste0(sim_dir, filename)
  )
  return(0)
}

### Simulation Case 7

sim_error_case7 <- function(max_time, interval, amp_noise, shape_noise, time_change, time_trend, n, run = 1, print_plots = FALSE) {
  require(tidyverse)
  require(pme)
  require(princurve)
  source("code/functions/sim_data.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  source("code/prinSurf_v3.R")

  time_vals <- seq(0, max_time, interval)
  sim_list <- lapply(
    time_vals,
    sim_data,
    case = 7,
    noise = 0.15,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    N = n,
    time_change = time_change,
    time_trend = time_trend
  )

  sim_df <- matrix(ncol = ncol(sim_list[[1]][[1]]))
  true_vals <- matrix(ncol = ncol(sim_list[[1]][[2]]))
  for (i in 1:length(sim_list)) {
    sim_df <- rbind(sim_df, sim_list[[i]][[1]])
    true_vals <- rbind(true_vals, sim_list[[i]][[2]])
  }
  sim_df <- sim_df[-1, ]
  true_vals <- true_vals[-1, ]


  lpme_result <- lpme(
    sim_df,
    2,
    smoothing_method = "spline",
    print_plots = print_plots,
    verbose = TRUE
  )
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)

  # lpme_result_gp <- lpme(
  #   sim_df,
  #   2,
  #   smoothing_method = "gp",
  #   print_plots = print_plots,
  #   verbose = TRUE
  # )
  # lpme_vals_gp <- calc_lpme_est(lpme_result_gp, sim_df)

  pme_result <- list()
  pme_vals <- list()
  principal_curve_result <- list()
  principal_curve_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 2, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
    principal_surface <- prinSurf(temp_data)
    surface_mse <- map(
      1:length(principal_surface),
      ~ principal_surface[[.x]]$MSE
    ) %>%
      unlist()
    opt_surface <- which.min(surface_mse)
    principal_curve_result[[t]] <- principal_surface[[opt_surface]]
    principal_curve_vals[[t]] <- cbind(time_vals[t], principal_surface[[opt_surface]]$PS)
  }
  pme_vals <- reduce(pme_vals, rbind)
  principal_curve_vals <- reduce(principal_curve_vals, rbind)

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], pme_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], lpme_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  # lpme_error_gp <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, ], lpme_vals_gp[.x, ])^2
  # ) %>%
  #   unlist() %>%
  #   mean()

  principal_curve_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], principal_curve_vals[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  p <- plot_ly(
    x = lpme_vals[, 2],
    y = lpme_vals[, 3],
    z = lpme_vals[, 4],
    frame = lpme_vals[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  ) %>%
    add_markers(
      x = pme_vals[, 2],
      y = pme_vals[, 3],
      z = pme_vals[, 4],
      frame = pme_vals[, 1]
    ) %>%
    add_markers(
      x = principal_curve_vals[, 2],
      y = principal_curve_vals[, 3],
      z = principal_curve_vals[, 4],
      frame = principal_curve_vals[, 1]
    ) %>%
    add_markers(
      x = true_vals[, 2],
      y = true_vals[, 3],
      z = true_vals[, 4],
      frame = true_vals[, 1]
    )

  sim_case7 <- list(
    df = sim_df,
    times = time_vals,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    n = n,
    lpme_result = lpme_result,
    # lpme_gp_result = lpme_result_gp,
    pme_results = pme_result,
    principal_curve_results = principal_curve_result,
    lpme_error = lpme_error,
    # lpme_gp_error = lpme_error_gp,
    pme_error = pme_error,
    principal_curve_error = principal_curve_error,
    plot = p
  )
  sim_dir <- "simulations/case7/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, recursive = TRUE)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_ampnoise_",
    str_pad(as.character(100 * amp_noise), 3, side = "left", pad = "0"),
    "_pernoise_",
    str_pad(as.character(100 * shape_noise), 3, side = "left", pad = "0"),
    "_n_",
    str_pad(as.character(n), 4, side = "left", pad = "0"),
    "_",
    time_trend,
    str_pad(as.character(100 * time_change), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case7,
    paste0(sim_dir, filename)
  )
  return(0)
}

### Simulation Case 8

sim_error_case8 <- function(max_time, interval, amp_noise, shape_noise, time_change, time_trend, n, run = 1, print_plots = FALSE) {
  require(tidyverse)
  require(pme)
  require(princurve)
  source("code/functions/sim_data.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  source("code/prinSurf_v3.R")


  time_vals <- seq(0, max_time, interval)
  sim_list <- lapply(
    time_vals,
    sim_data,
    case = 8,
    noise = 0.15,
    amp_noise = amp_noise / 25,
    period_noise = shape_noise / 25,
    N = n,
    time_change = time_change,
    time_trend = time_trend
  )

  sim_df <- matrix(ncol = ncol(sim_list[[1]][[1]]))
  true_vals <- matrix(ncol = ncol(sim_list[[1]][[2]]))
  for (i in 1:length(sim_list)) {
    sim_df <- rbind(sim_df, sim_list[[i]][[1]])
    true_vals <- rbind(true_vals, sim_list[[i]][[2]])
  }
  sim_df <- sim_df[-1, ]
  true_vals <- true_vals[-1, ]

  pol <- sim_df[, -c(1, 4)] %>%
    cart2pol()
  sim_df <- cbind(sim_df, pol)
  sim_df <- sim_df[, -5]

  lpme_result <- lpme(
    sim_df,
    2,
    smoothing_method = "spline",
    print_plots = print_plots,
    verbose = TRUE
  )
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)

  # lpme_result_gp <- lpme(
  #   sim_df,
  #   2,
  #   smoothing_method = "gp",
  #   print_plots = print_plots,
  #   verbose = TRUE
  # )
  # lpme_vals_gp <- calc_lpme_est(lpme_result_gp, sim_df)

  pme_result <- list()
  pme_vals <- list()
  principal_curve_result <- list()
  principal_curve_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 2, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
    principal_surface <- prinSurf(temp_data)
    surface_mse <- map(
      1:length(principal_surface),
      ~ principal_surface[[.x]]$MSE
    ) %>%
      unlist()
    opt_surface <- which.min(surface_mse)
    principal_curve_result[[t]] <- principal_surface[[opt_surface]]
    principal_curve_vals[[t]] <- cbind(time_vals[t], principal_surface[[opt_surface]]$PS)
  }
  pme_vals <- reduce(pme_vals, rbind)
  principal_curve_vals <- reduce(principal_curve_vals, rbind)

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], pme_vals[.x, 1:4])^2
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], lpme_vals[.x, 1:4])^2
  ) %>%
    unlist() %>%
    mean()

  # lpme_error_gp <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, ], lpme_vals_gp[.x, 1:4])^2
  # ) %>%
  #   unlist() %>%
  #   mean()

  principal_curve_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], principal_curve_vals[.x, 1:4])^2
  ) %>%
    unlist() %>%
    mean()

  p <- plot_ly(
    x = lpme_vals[, 2],
    y = lpme_vals[, 3],
    z = lpme_vals[, 4],
    frame = lpme_vals[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  ) %>%
    add_markers(
      x = pme_vals[, 2],
      y = pme_vals[, 3],
      z = pme_vals[, 4],
      frame = pme_vals[, 1]
    ) %>%
    add_markers(
      x = principal_curve_vals[, 2],
      y = principal_curve_vals[, 3],
      z = principal_curve_vals[, 4],
      frame = principal_curve_vals[, 1]
    ) %>%
    add_markers(
      x = true_vals[, 2],
      y = true_vals[, 3],
      z = true_vals[, 4],
      frame = true_vals[, 1]
    )


  sim_case8 <- list(
    df = sim_df,
    times = time_vals,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    n = n,
    lpme_result = lpme_result,
    # lpme_gp_result = lpme_result_gp,
    pme_results = pme_result,
    principal_curve_results = principal_curve_result,
    lpme_error = lpme_error,
    # lpme_gp_error = lpme_error_gp,
    pme_error = pme_error,
    principal_curve_error = principal_curve_error,
    plot = p
  )
  sim_dir <- "simulations/case8/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, recursive = TRUE)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_ampnoise_",
    str_pad(as.character(100 * amp_noise), 3, side = "left", pad = "0"),
    "_pernoise_",
    str_pad(as.character(100 * shape_noise), 3, side = "left", pad = "0"),
    "_n_",
    str_pad(as.character(n), 4, side = "left", pad = "0"),
    "_",
    time_trend,
    str_pad(as.character(100 * time_change), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case8,
    paste0(sim_dir, filename)
  )
  return(0)
}

### Simulation Case 9

sim_error_case9 <- function(max_time, interval, amp_noise, shape_noise, time_change, time_trend, n, run = 1, print_plots = FALSE) {
  require(tidyverse)
  require(pme)
  require(princurve)
  source("code/functions/sim_data.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  source("code/prinSurf_v3.R")


  time_vals <- seq(0, max_time, interval)
  sim_list <- lapply(
    time_vals,
    sim_data,
    case = 9,
    noise = 0.15,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    N = n,
    time_change = time_change,
    time_trend = time_trend
  )

  sim_df <- matrix(ncol = ncol(sim_list[[1]][[1]]))
  true_vals <- matrix(ncol = ncol(sim_list[[1]][[2]]))
  for (i in 1:length(sim_list)) {
    sim_df <- rbind(sim_df, sim_list[[i]][[1]])
    true_vals <- rbind(true_vals, sim_list[[i]][[2]])
  }
  sim_df <- sim_df[-1, ]
  true_vals <- true_vals[-1, ]

  sph <- sim_df[, -1] %>%
    cart2sph()
  sim_df <- cbind(sim_df, sph)


  lpme_result <- lpme(
    sim_df,
    2,
    smoothing_method = "spline",
    print_plots = print_plots,
    verbose = TRUE
  )
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)

  # lpme_result_gp <- lpme(
  #   sim_df,
  #   2,
  #   smoothing_method = "gp",
  #   print_plots = print_plots,
  #   verbose = TRUE
  # )
  # lpme_vals_gp <- calc_lpme_est(lpme_result_gp, sim_df)

  pme_result <- list()
  pme_vals <- list()
  principal_curve_result <- list()
  principal_curve_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 2, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
    principal_surface <- prinSurf(temp_data)
    surface_mse <- map(
      1:length(principal_surface),
      ~ principal_surface[[.x]]$MSE
    ) %>%
      unlist()
    opt_surface <- which.min(surface_mse)
    principal_curve_result[[t]] <- principal_surface[[opt_surface]]
    principal_curve_vals[[t]] <- cbind(time_vals[t], principal_surface[[opt_surface]]$PS)
  }
  pme_vals <- reduce(pme_vals, rbind)
  principal_curve_vals <- reduce(principal_curve_vals, rbind)

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], pme_vals[.x, 1:4])^2
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], lpme_vals[.x, 1:4])^2
  ) %>%
    unlist() %>%
    mean()

  # lpme_error_gp <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, ], lpme_vals_gp[.x, 1:4])^2
  # ) %>%
  #   unlist() %>%
  #   mean()

  principal_curve_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], principal_curve_vals[.x, 1:4])^2
  ) %>%
    unlist() %>%
    mean()

  p <- plot_ly(
    x = lpme_vals[, 2],
    y = lpme_vals[, 3],
    z = lpme_vals[, 4],
    frame = lpme_vals[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  ) %>%
    add_markers(
      x = pme_vals[, 2],
      y = pme_vals[, 3],
      z = pme_vals[, 4],
      frame = pme_vals[, 1]
    ) %>%
    add_markers(
      x = principal_curve_vals[, 2],
      y = principal_curve_vals[, 3],
      z = principal_curve_vals[, 4],
      frame = principal_curve_vals[, 1]
    ) %>%
    add_markers(
      x = true_vals[, 2],
      y = true_vals[, 3],
      z = true_vals[, 4],
      frame = true_vals[, 1]
    )

  sim_case9 <- list(
    df = sim_df,
    times = time_vals,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    n = n,
    lpme_result = lpme_result,
    # lpme_gp_result = lpme_result_gp,
    pme_results = pme_result,
    principal_curve_results = principal_curve_result,
    lpme_error = lpme_error,
    # lpme_gp_error = lpme_error_gp,
    pme_error = pme_error,
    principal_curve_error = principal_curve_error,
    plot = p
  )
  sim_dir <- "simulations/case9/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, recursive = TRUE)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_ampnoise_",
    str_pad(as.character(100 * amp_noise), 3, side = "left", pad = "0"),
    "_pernoise_",
    str_pad(as.character(100 * shape_noise), 3, side = "left", pad = "0"),
    "_n_",
    str_pad(as.character(n), 4, side = "left", pad = "0"),
    "_",
    time_trend,
    str_pad(as.character(100 * time_change), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case9,
    paste0(sim_dir, filename)
  )
  return(0)
}

### Simulation Case 10

sim_error_case10 <- function(max_time, interval, amp_noise, shape_noise, time_change, time_trend, n, run = 1, print_plots = FALSE) {
  require(tidyverse)
  require(pme)
  require(princurve)
  source("code/functions/sim_data.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  source("code/prinSurf_v3.R")


  time_vals <- seq(0, max_time, interval)
  sim_list <- lapply(
    time_vals,
    sim_data,
    case = 10,
    noise = 0.15,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    N = n,
    time_change = time_change,
    time_trend = time_trend
  )

  sim_df <- matrix(ncol = ncol(sim_list[[1]][[1]]))
  true_vals <- matrix(ncol = ncol(sim_list[[1]][[2]]))
  for (i in 1:length(sim_list)) {
    sim_df <- rbind(sim_df, sim_list[[i]][[1]])
    true_vals <- rbind(true_vals, sim_list[[i]][[2]])
  }
  sim_df <- sim_df[-1, ]
  true_vals <- true_vals[-1, ]

  sph <- sim_df[, -1] %>%
    cart2sph()
  sim_df <- cbind(sim_df, sph)


  lpme_result <- lpme(
    sim_df,
    2,
    smoothing_method = "spline",
    print_plots = print_plots,
    verbose = TRUE
  )
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)

  # lpme_result_gp <- lpme(
  #   sim_df,
  #   2,
  #   smoothing_method = "gp",
  #   print_plots = print_plots,
  #   verbose = TRUE
  # )
  # lpme_vals_gp <- calc_lpme_est(lpme_result_gp, sim_df)

  pme_result <- list()
  pme_vals <- list()
  principal_curve_result <- list()
  principal_curve_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 2, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
    principal_surface <- prinSurf(temp_data)
    surface_mse <- map(
      1:length(principal_surface),
      ~ principal_surface[[.x]]$MSE
    ) %>%
      unlist()
    opt_surface <- which.min(surface_mse)
    principal_curve_result[[t]] <- principal_surface[[opt_surface]]
    principal_curve_vals[[t]] <- cbind(time_vals[t], principal_surface[[opt_surface]]$PS)
  }
  pme_vals <- reduce(pme_vals, rbind)
  principal_curve_vals <- reduce(principal_curve_vals, rbind)

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], pme_vals[.x, 1:4])^2
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], lpme_vals[.x, 1:4])^2
  ) %>%
    unlist() %>%
    mean()

  # lpme_error_gp <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, ], lpme_vals_gp[.x, 1:4])^2
  # ) %>%
  #   unlist() %>%
  #   mean()

  principal_curve_error <- map(
    1:nrow(true_vals),
    ~ dist_euclidean(true_vals[.x, ], principal_curve_vals[.x, 1:4])^2
  ) %>%
    unlist() %>%
    mean()

  p <- plot_ly(
    x = lpme_vals[, 2],
    y = lpme_vals[, 3],
    z = lpme_vals[, 4],
    frame = lpme_vals[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  ) %>%
    add_markers(
      x = pme_vals[, 2],
      y = pme_vals[, 3],
      z = pme_vals[, 4],
      frame = pme_vals[, 1]
    ) %>%
    add_markers(
      x = principal_curve_vals[, 2],
      y = principal_curve_vals[, 3],
      z = principal_curve_vals[, 4],
      frame = principal_curve_vals[, 1]
    ) %>%
    add_markers(
      x = true_vals[, 2],
      y = true_vals[, 3],
      z = true_vals[, 4],
      frame = true_vals[, 1]
    )


  sim_case10 <- list(
    df = sim_df,
    times = time_vals,
    amp_noise = amp_noise,
    period_noise = shape_noise,
    n = n,
    lpme_result = lpme_result,
    # lpme_gp_result = lpme_result_gp,
    pme_results = pme_result,
    principal_curve_results = principal_curve_result,
    lpme_error = lpme_error,
    # lpme_gp_error = lpme_error_gp,
    pme_error = pme_error,
    principal_curve_error = principal_curve_error,
    plot = p
  )
  sim_dir <- "simulations/case10/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, recursive = TRUE)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_ampnoise_",
    str_pad(as.character(100 * amp_noise), 3, side = "left", pad = "0"),
    "_pernoise_",
    str_pad(as.character(100 * shape_noise), 3, side = "left", pad = "0"),
    "_n_",
    str_pad(as.character(n), 4, side = "left", pad = "0"),
    "_",
    time_trend,
    str_pad(as.character(100 * time_change), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case10,
    paste0(sim_dir, filename)
  )
  return(0)
}

amp_noise_vals <- seq(0, 2, 0.25)
shape_noise_vals <- seq(0, 2, 0.25)
max_times <- c(1, 2, 5, 8, 10)
intervals <- c(0.1, 0.25, 0.5, 1)
n_vals <- 10^(3:4)
replicates <- 1:4
case <- 1:10
time_changes <- seq(0, 2, 0.25)
time_trends <- c("constant", "linear", "quadratic", "sinusoidal")

param_grid <- expand.grid(
  amp_noise_vals,
  shape_noise_vals,
  n_vals,
  intervals,
  max_times,
  replicates,
  case,
  time_changes,
  time_trends
)

param_grid <- param_grid[param_grid[, 7] != 4, ]

# plan(multisession, workers = availableCores() - 2)
# plan(multicore, workers = availableCores() - 2)
plan(sequential)
set.seed(21986)
# pb <- progress_bar$new(total = nrow(param_grid))
# errors <- map(
#   sample(1:nrow(param_grid), nrow(param_grid)),
#   # 1:nrow(param_grid),
#   ~ {
#       if (param_grid[.x, 7] == 1) {
#         try(
#           sim_error_case1(
#             param_grid[.x, 5],
#             param_grid[.x, 4],
#             param_grid[.x, 1],
#             param_grid[.x, 2],
#             param_grid[.x, 8],
#             param_grid[.x, 9],
#             param_grid[.x, 3],
#             param_grid[.x, 6],
#             print_plots = TRUE
#           )
#         )
#       } else if (param_grid[.x, 7] == 2) {
#         try(
#           sim_error_case2(
#             param_grid[.x, 5],
#             param_grid[.x, 4],
#             param_grid[.x, 1],
#             param_grid[.x, 2],
#             param_grid[.x, 8],
#             param_grid[.x, 9],
#             param_grid[.x, 3],
#             param_grid[.x, 6],
#             print_plots = TRUE
#           )
#         )
#       } else if (param_grid[.x, 7] == 3) {
#         try(
#           sim_error_case3(
#             param_grid[.x, 5],
#             param_grid[.x, 4],
#             param_grid[.x, 1],
#             param_grid[.x, 2],
#             param_grid[.x, 8],
#             param_grid[.x, 9],
#             param_grid[.x, 3],
#             param_grid[.x, 6],
#             print_plots = TRUE
#           )
#         )
#       } else if (param_grid[.x, 7] == 5) {
#         try(
#           sim_error_case5(
#             param_grid[.x, 5],
#             param_grid[.x, 4],
#             param_grid[.x, 1],
#             param_grid[.x, 2],
#             param_grid[.x, 8],
#             param_grid[.x, 9],
#             param_grid[.x, 3],
#             param_grid[.x, 6],
#             print_plots = TRUE
#           )
#         )
#       } else if (param_grid[.x, 7] == 6) {
#         try(
#           sim_error_case6(
#             param_grid[.x, 5],
#             param_grid[.x, 4],
#             param_grid[.x, 1],
#             param_grid[.x, 2],
#             param_grid[.x, 8],
#             param_grid[.x, 9],
#             param_grid[.x, 3],
#             param_grid[.x, 6],
#             print_plots = TRUE
#           )
#         )
#       } else if (param_grid[.x, 7] == 7) {
#         try(
#           sim_error_case7(
#             param_grid[.x, 5],
#             param_grid[.x, 4],
#             param_grid[.x, 1],
#             param_grid[.x, 2],
#             param_grid[.x, 8],
#             param_grid[.x, 9],
#             param_grid[.x, 3],
#             param_grid[.x, 6],
#             print_plots = TRUE
#           )
#         )
#       } else if (param_grid[.x, 7] == 8) {
#         try(
#           sim_error_case8(
#             param_grid[.x, 5],
#             param_grid[.x, 4],
#             param_grid[.x, 1],
#             param_grid[.x, 2],
#             param_grid[.x, 8],
#             param_grid[.x, 9],
#             param_grid[.x, 3],
#             param_grid[.x, 6],
#             print_plots = TRUE
#           )
#         )
#       } else if (param_grid[.x, 7] == 9) {
#         try(
#           sim_error_case9(
#             param_grid[.x, 5],
#             param_grid[.x, 4],
#             param_grid[.x, 1],
#             param_grid[.x, 2],
#             param_grid[.x, 8],
#             param_grid[.x, 9],
#             param_grid[.x, 3],
#             param_grid[.x, 6],
#             print_plots = TRUE
#           )
#         )
#       } else if (param_grid[.x, 7] == 10) {
#         try(
#           sim_error_case10(
#             param_grid[.x, 5],
#             param_grid[.x, 4],
#             param_grid[.x, 1],
#             param_grid[.x, 2],
#             param_grid[.x, 8],
#             param_grid[.x, 9],
#             param_grid[.x, 3],
#             param_grid[.x, 6],
#             print_plots = TRUE
#           )
#         )
#       }
#     }
# )

errors <- future_map(
  sample(1:nrow(param_grid), nrow(param_grid)),
  ~ {
      if (param_grid[.x, 7] == 1) {
        try(
          sim_error_case1(
            param_grid[.x, 5],
            param_grid[.x, 4],
            param_grid[.x, 1],
            param_grid[.x, 2],
            param_grid[.x, 8],
            param_grid[.x, 9],
            param_grid[.x, 3],
            param_grid[.x, 6],
            print_plots = FALSE
          )
        )
      } else if (param_grid[.x, 7] == 2) {
        try(
          sim_error_case2(
            param_grid[.x, 5],
            param_grid[.x, 4],
            param_grid[.x, 1],
            param_grid[.x, 2],
            param_grid[.x, 8],
            param_grid[.x, 9],
            param_grid[.x, 3],
            param_grid[.x, 6],
            print_plots = FALSE
          )
        )
      } else if (param_grid[.x, 7] == 3) {
        try(
          sim_error_case3(
            param_grid[.x, 5],
            param_grid[.x, 4],
            param_grid[.x, 1],
            param_grid[.x, 2],
            param_grid[.x, 8],
            param_grid[.x, 9],
            param_grid[.x, 3],
            param_grid[.x, 6],
            print_plots = FALSE
          )
        )
      } else if (param_grid[.x, 7] == 5) {
        try(
          sim_error_case5(
            param_grid[.x, 5],
            param_grid[.x, 4],
            param_grid[.x, 1],
            param_grid[.x, 2],
            param_grid[.x, 8],
            param_grid[.x, 9],
            param_grid[.x, 3],
            param_grid[.x, 6],
            print_plots = FALSE
          )
        )
      } else if (param_grid[.x, 7] == 6) {
        try(
          sim_error_case6(
            param_grid[.x, 5],
            param_grid[.x, 4],
            param_grid[.x, 1],
            param_grid[.x, 2],
            param_grid[.x, 8],
            param_grid[.x, 9],
            param_grid[.x, 3],
            param_grid[.x, 6],
            print_plots = FALSE
          )
        )
      } else if (param_grid[.x, 7] == 7) {
        try(
          sim_error_case7(
            param_grid[.x, 5],
            param_grid[.x, 4],
            param_grid[.x, 1],
            param_grid[.x, 2],
            param_grid[.x, 8],
            param_grid[.x, 9],
            param_grid[.x, 3],
            param_grid[.x, 6],
            print_plots = FALSE
          )
        )
      } else if (param_grid[.x, 7] == 8) {
        try(
          sim_error_case8(
            param_grid[.x, 5],
            param_grid[.x, 4],
            param_grid[.x, 1],
            param_grid[.x, 2],
            param_grid[.x, 8],
            param_grid[.x, 9],
            param_grid[.x, 3],
            param_grid[.x, 6],
            print_plots = FALSE
          )
        )
      } else if (param_grid[.x, 7] == 9) {
        try(
          sim_error_case9(
            param_grid[.x, 5],
            param_grid[.x, 4],
            param_grid[.x, 1],
            param_grid[.x, 2],
            param_grid[.x, 8],
            param_grid[.x, 9],
            param_grid[.x, 3],
            param_grid[.x, 6],
            print_plots = FALSE
          )
        )
      } else if (param_grid[.x, 7] == 10) {
        try(
          sim_error_case10(
            param_grid[.x, 5],
            param_grid[.x, 4],
            param_grid[.x, 1],
            param_grid[.x, 2],
            param_grid[.x, 8],
            param_grid[.x, 9],
            param_grid[.x, 3],
            param_grid[.x, 6],
            print_plots = FALSE
          )
        )
      }
    },
  .progress = TRUE,
  .options = furrr_options(seed = TRUE, scheduling = Inf)
)
