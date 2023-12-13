library(tidyverse)
library(plotly)
library(pracma)
library(foreach)
library(doParallel)
library(furrr)
library(progress)
library(pme)
library(princurve)
library(Rfast)

source("code/functions/sim_data.R")
source("code/functions/calc_pme_est.R")
source("code/functions/calc_lpme_est.R")
source("code/prinSurf_v3.R")
source("code/functions/estimate_volume.R")

est_sphere_vol <- function(max_time, interval, amp_noise, r, pct_missingness, run) {

  time_vals <- seq(0, max_time, interval)
  sim_list <- lapply(
    time_vals,
    sim_data,
    case = 9,
    noise = 0.05,
    amp_noise = amp_noise,
    period_noise = 0,
    N = 1000,
    time_change = 0,
    time_trend = "constant"
  )

  sim_df <- matrix(ncol = ncol(sim_list[[1]][[1]]))
  true_vals <- matrix(ncol = ncol(sim_list[[1]][[2]]))
  for (i in 1:length(sim_list)) {
    sim_df <- rbind(sim_df, sim_list[[i]][[1]])
    true_vals <- rbind(true_vals, sim_list[[i]][[2]])
  }
  sim_df <- sim_df[-1, ]
  true_vals <- true_vals[-1, ]

  for (col_idx in 2:ncol(sim_df)) {
    sim_df[, col_idx] <- r * sim_df[, col_idx]
    true_vals[, col_idx] <- r * true_vals[, col_idx]
  }

  true_sph <- cart2sph(true_vals[, -1])

  true_r <- vector()
  true_volume <- vector()
  for (time_idx in 1:length(unique(true_vals[, 1]))) {
    time_val <- unique(true_vals[, 1])[time_idx]
    true_r[time_idx] <- mean(true_sph[true_vals[, 1] == time_val, 3])
    true_volume[time_idx] <- (4 / 3) * pi * true_r[time_idx]^3
  }

  col_maxs <- vector()
  for (col_idx in 1:ncol(sim_df[, -(5:7)])) {
    col_max <- max(abs(sim_df[, col_idx]))
    col_maxs[col_idx] <- col_max
    sim_df[, col_idx] <- sim_df[, col_idx] / col_max
    true_vals[, col_idx] <- true_vals[, col_idx] / col_max
  }

  sph <- sim_df[, -1] %>%
    cart2sph()
  true_sph <- true_vals[, -1] %>%
    cart2sph()
  sim_df <- cbind(sim_df, sph)
  true_vals <- cbind(true_vals, true_sph)

  time_vals <- unique(sim_df[, 1])

  lpme_isomap <- lpme(
    sim_df,
    2,
    smoothing_method = "spline",
    print_plots = FALSE,
    verbose = TRUE,
    init = "first",
    initialization_algorithm = "isomap"
  )

  lpme_dm <- lpme(
    sim_df,
    2,
    smoothing_method = "spline",
    print_plots = FALSE,
    verbose = TRUE,
    init = "first",
    initialization_algorithm = "diffusion_maps"
  )

  lpme_le <- lpme(
    sim_df,
    2,
    smoothing_method = "spline",
    print_plots = FALSE,
    verbose = TRUE,
    init = "first",
    initialization_algorithm = "laplacian_eigenmaps"
  )

  lpme_vals_isomap <- calc_lpme_est(lpme_isomap, sim_df)
  lpme_vals_dm <- calc_lpme_est(lpme_dm, sim_df)
  lpme_vals_le <- calc_lpme_est(lpme_le, sim_df)

  pme_isomap <- list()
  pme_dm <- list()
  pme_le <- list()
  pme_vals_isomap <- list()
  pme_vals_dm <- list()
  pme_vals_le <- list()

  principal_curve_result <- list()
  principal_curve_vals <- list()

  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_isomap[[t]] <- pme(
      temp_data,
      d = 2,
      initialization_algorithm = "isomap",
      initialization_type = "centers"
    )
    pme_dm[[t]] <- pme(
      temp_data,
      d = 2,
      initialization_algorithm = "diffusion_maps",
      initialization_type = "centers"
    )
    pme_le[[t]] <- pme(
      temp_data,
      d = 2,
      initialization_algorithm = "laplacian_eigenmaps",
      initialization_type = "centers"
    )
    pme_vals_isomap[[t]] <- cbind(
      time_vals[t],
      calc_pme_est(pme_isomap[[t]], temp_data)
    )
    pme_vals_dm[[t]] <- cbind(
      time_vals[t],
      calc_pme_est(pme_dm[[t]], temp_data)
    )
    pme_vals_le[[t]] <- cbind(
      time_vals[t],
      calc_pme_est(pme_le[[t]], temp_data)
    )

    # principal_surface <- prinSurf(temp_data)
    # surface_mse <- map(
    #   1:length(principal_surface),
    #   ~ principal_surface[[.x]]$MSE
    # ) %>%
    #   unlist()
    # opt_surface <- which.min(surface_mse)
    # principal_curve_result[[t]] <- principal_surface[[opt_surface + 2]]
    # principal_curve_vals[[t]] <- cbind(time_vals[t], principal_surface[[opt_surface + 2]]$PS)
  }
  pme_vals_isomap <- reduce(pme_vals_isomap, rbind)
  pme_vals_dm <- reduce(pme_vals_dm, rbind)
  pme_vals_le <- reduce(pme_vals_le, rbind)
  # principal_curve_vals <- reduce(principal_curve_vals, rbind)

  # pme_error <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, 1:4], pme_vals[.x, 1:4])^2
  # ) %>%
  #   unlist() %>%
  #   mean()
  #
  # lpme_error <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, 1:4], lpme_vals[.x, 1:4])^2
  # ) %>%
  #   unlist() %>%
  #   mean()
  #
  # data_error <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, 1:4], sim_df[.x, 1:4])^2
  # ) %>%
  #   unlist() %>%
  #   mean()

  # lpme_error_gp <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, ], lpme_vals_gp[.x, 1:4])^2
  # ) %>%
  #   unlist() %>%
  #   mean()

  # principal_curve_error <- map(
  #   1:nrow(true_vals),
  #   ~ dist_euclidean(true_vals[.x, 1:4], principal_curve_vals[.x, 1:4])^2
  # ) %>%
  #   unlist() %>%
  #   mean()


  sim_data_rescaled <- sim_df
  lpme_isomap_rescaled <- lpme_vals_isomap
  lpme_dm_rescaled <- lpme_vals_dm
  lpme_le_rescaled <- lpme_vals_le
  pme_isomap_rescaled <- pme_vals_isomap
  pme_dm_rescaled <- pme_vals_dm
  pme_le_rescaled <- pme_vals_le

  for (col_idx in 2:ncol(sim_df)) {
    sim_data_rescaled[, col_idx] <- round(sim_data_rescaled[, col_idx] * col_maxs[col_idx])
    lpme_isomap_rescaled[, col_idx] <- round(lpme_isomap_rescaled[, col_idx] * col_maxs[col_idx])
    pme_isomap_rescaled[, col_idx] <- round(pme_isomap_rescaled[, col_idx] * col_maxs[col_idx])
  }
  sim_data_rescaled <- sim_data_rescaled[, -(5:7)]
  lpme_isomap_rescaled <- lpme_isomap_rescaled[, -(5:7)]
  pme_isomap_rescaled <- pme_isomap_rescaled[, -(5:7)]

  est_volume_data <- list()
  est_volume_lpme_isomap <- list()
  est_volume_lpme_dm <- list()
  est_volume_lpme_le <- list()
  est_volume_pme_isomap <- list()
  est_volume_pme_dm <- list()
  est_volume_pme_le <- list()
  for (time_idx in 1:length(unique(sim_data_rescaled[, 1]))) {
    temp_data <- sim_data_rescaled[sim_data_rescaled[, 1] == time_vals[time_idx], -1]
    missing_vals <- sample(1:nrow(temp_data), round(pct_missingness * nrow(temp_data)))
    if (length(missing_vals) > 0) {
      temp_data <- temp_data[-missing_vals, ]
    }
    est_volume_data[[time_idx]] <- estimate_volume(temp_data, 1)

    temp_lpme_isomap <- lpme_isomap_rescaled[lpme_isomap_rescaled[, 1] == time_vals[time_idx], -1]
    est_volume_lpme_isomap[[time_idx]] <- estimate_volume(temp_lpme_isomap, 1)

    temp_lpme_dm <- lpme_dm_rescaled[lpme_dm_rescaled[, 1] == time_vals[time_idx], -1]
    est_volume_lpme_dm[[time_idx]] <- estimate_volume(temp_lpme_dm, 1)

    temp_lpme_le <- lpme_le_rescaled[lpme_le_rescaled[, 1] == time_vals[time_idx], -1]
    est_volume_lpme_le[[time_idx]] <- estimate_volume(temp_lpme_le, 1)

    temp_pme_isomap <- pme_isomap_rescaled[pme_isomap_rescaled[, 1] == time_vals[time_idx], -1]
    est_volume_pme_isomap[[time_idx]] <- estimate_volume(temp_pme_isomap, 1)

    temp_pme_dm <- pme_dm_rescaled[pme_dm_rescaled[, 1] == time_vals[time_idx], -1]
    est_volume_pme_dm[[time_idx]] <- estimate_volume(temp_pme_dm, 1)

    temp_pme_le <- pme_le_rescaled[pme_le_rescaled[, 1] == time_vals[time_idx], -1]
    est_volume_pme_le[[time_idx]] <- estimate_volume(temp_pme_le, 1)
  }

  est_volume_points_data <- map(est_volume_data, ~ .x[[2]])
  est_volume_data <- map(est_volume_data, ~ .x[[1]]) %>%
    unlist()
  est_volume_points_lpme_isomap <- map(est_volume_lpme_isomap, ~ .x[[2]])
  est_volume_lpme_isomap <- map(est_volume_lpme_isomap, ~ .x[[1]]) %>%
    unlist()
  est_volume_points_lpme_dm <- map(est_volume_lpme_dm, ~ .x[[2]])
  est_volume_lpme_dm <- map(est_volume_lpme_dm, ~ .x[[1]]) %>%
    unlist()
  est_volume_points_lpme_le <- map(est_volume_lpme_le, ~ .x[[2]])
  est_volume_lpme_le <- map(est_volume_lpme_le, ~ .x[[1]]) %>%
    unlist()
  est_volume_points_pme_isomap <- map(est_volume_pme_isomap, ~ .x[[2]])
  est_volume_pme_isomap <- map(est_volume_pme_isomap, ~ .x[[1]]) %>%
    unlist()
  est_volume_points_pme_dm <- map(est_volume_pme_dm, ~ .x[[2]])
  est_volume_pme_dm <- map(est_volume_pme_dm, ~ .x[[1]]) %>%
    unlist()
  est_volume_points_pme_le <- map(est_volume_pme_le, ~ .x[[2]])
  est_volume_pme_le <- map(est_volume_pme_le, ~ .x[[1]]) %>%
    unlist()

  volume_df <- tibble(
    time = time_vals,
    amp_noise = amp_noise,
    r = r,
    pct_missing = pct_missingness,
    true_volume = true_volume,
    data_volume = est_volume_data,
    lpme_isomap_volume = est_volume_lpme_isomap,
    lpme_dm_volume = est_volume_lpme_dm,
    lpme_le_volume = est_volume_lpme_le,
    pme_isomap_volume = est_volume_pme_isomap,
    pme_dm_volume = est_volume_pme_dm,
    pme_le_volume = est_volume_pme_le
  ) %>%
    mutate(
      data_error = (data_volume - true_volume) / true_volume,
      lpme_isomap_error = (lpme_isomap_volume - true_volume) / true_volume,
      lpme_dm_error = (lpme_dm_volume - true_volume) / true_volume,
      lpme_le_error = (lpme_le_volume - true_volume) / true_volume,
      pme_isomap_error = (pme_isomap_volume - true_volume) / true_volume,
      pme_dm_error = (pme_dm_volume - true_volume) / true_volume,
      pme_le_error = (pme_le_volume - true_volume) / true_volume
    )
  volume_out <- list(
    volume_data = volume_df,
    data_error_mean = mean(volume_df$data_error),
    data_error_sd = sd(volume_df$data_error),
    lpme_isomap_error_mean = mean(volume_df$lpme_isomap_error),
    lpme_isomap_error_sd = sd(volume_df$lpme_isomap_error),
    lpme_dm_error_mean = mean(volume_df$lpme_dm_error),
    lpme_dm_error_sd = sd(volume_df$lpme_dm_error),
    lpme_le_error_mean = mean(volume_df$lpme_le_error),
    lpme_le_error_sd = sd(volume_df$lpme_le_error),
    pme_isomap_error_mean = mean(volume_df$pme_isomap_error),
    pme_isomap_error_sd = sd(volume_df$pme_isomap_error),
    pme_dm_error_mean = mean(volume_df$pme_dm_error),
    pme_dm_error_sd = sd(volume_df$pme_dm_error),
    pme_le_error_mean = mean(volume_df$pme_le_error),
    pme_le_error_sd = sd(volume_df$pme_le_error),
    duration = max_time,
    interval = interval,
    noise = amp_noise,
    r = r,
    pct_missing = pct_missingness,
    true_volume = true_volume
  )
  volume_dir <- "simulations/volume"
  if (!dir.exists(volume_dir)) {
    dir.create(volume_dir)
  }

  saveRDS(
    volume_out,
    paste0(
      volume_dir,
      "/duration_",
      str_pad(max_time, 2, pad = "0"),
      "_interval_",
      str_pad(round(interval * 100), 2, pad = "0"),
      "_r_",
      str_pad(r, 2, pad = "0"),
      "_noise_",
      str_pad(round(amp_noise * 100), 2, pad = "0"),
      "_missing_",
      str_pad(round(pct_missingness * 100), 2, pad = "0"),
      "_run",
      str_pad(run, 2, pad = "0"),
      ".rds"
    )
  )
  return(0)
}

durations <- c(1, 2, 5)
intervals <- c(0.1, 0.25, 0.5)
noise_vals <- c(0, 0.1, 0.25, 0.5)
r_vals <- c(1, 5, 10, 25)
missingness <- c(0, 0.05, 0.1, 0.25)
runs <- 1:10

param_grid <- expand.grid(
  durations,
  intervals,
  noise_vals,
  r_vals,
  missingness,
  runs
)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

set.seed(94713)
volumes <- foreach(
  # x = sample(1:nrow(param_grid), nrow(param_grid)),
  duration = param_grid[, 1],
  interval = param_grid[, 2],
  noise = param_grid[, 3],
  r = param_grid[, 4],
  missing_pct = param_grid[, 5],
  run = param_grid[, 6],
  .inorder = FALSE,
  .export = c("sim_data", "calc_pme_est", "calc_lpme_est", "cart2sph", "estimate_volume"),
  .packages = c("tidyverse", "pme", "princurve", "plotly", "doParallel")
) %dopar% {
    est_sphere_vol(duration, interval, noise, r, missing_pct, run)
  }

stopCluster(cl)
