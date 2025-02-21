estimate_volume_case8 <- function(sim, n_points) {
  require(Rfast, warn.conflicts = FALSE, quietly = TRUE)

  source("functions/interior_identification.R")
  sim_data <- sim$processed_data$df
  time_points <- sim_data |>
    select(time) |>
    unlist() |>
    unique()

  observed_times <- sim$data$df |>
    select(time) |>
    unlist() |>
    unique() |>
    sort()

  sim_mat <- sim$processed_data$df_observed |>
    select(time, contains("X")) |>
    as.matrix()

  radii <- map(
    seq_along(sim$data$amplitude_values),
    ~ sim$data$amplitude_values[[.x]][1]
  ) |>
    reduce(c)

  observed_volumes <- (4 / 3) * pi * radii^3

  time_adjustments <- case_when(
    time_trend_val == "constant" ~ 0 * time_points,
    time_trend_val == "linear" ~ time_points,
    time_trend_val == "quadratic" ~ time_points^2
  )

# scale time adjustments to be proportion of 1
  if (max(time_adjustments) == 0) {
    time_adjustments_scaled <- time_adjustments * time_change_val * 0
  } else {
    time_adjustments_scaled <- (
      time_adjustments / max(time_adjustments)
    ) * time_change_val
  }
# subtract from initial amplitude to replicate effects of atrophy
  true_radii <- 1 - time_adjustments_scaled
  true_volumes <- (4 / 3) * pi * true_radii^3 

  ## Calcalate volumes from partitioned models

  lpme_volumes <- vector(mode = "numeric", length = length(time_points))
  pme_volumes <- vector(mode = "numeric", length = length(time_points))

  lpme_part_volumes <- vector(mode = "numeric", length = length(time_points))
  pme_part_volumes <- vector(mode = "numeric", length = length(time_points))

  scaling_factors <- sim$data$df |>
    select(contains("X")) |>
    select(!contains("true")) |>
    colMinsMaxs()
  x_scale <- max(abs(scaling_factors[, 1]))
  y_scale <- max(abs(scaling_factors[, 2]))
  z_scale <- max(abs(scaling_factors[, 3]))



  lpme_est <- sim$models$lpme$lpme
  lpme_params <- lpme_est$parameterization_list[[which.min(lpme_est$msd)]] 


  for (time_idx in seq_along(time_points)) {
    observed_data <- sim$data$df |>
      filter(time == observed_times[time_idx]) |>
      select(contains("X")) |>
      select(!contains("true")) |>
      as.matrix()
    obs_limits <- colMinsMaxs(observed_data)
    limit_scaler <- 0.05
    obs_x_range <- obs_limits[2, 1] - obs_limits[1, 1]
    obs_y_range <- obs_limits[2, 2] - obs_limits[1, 2]
    obs_z_range <- obs_limits[2, 3] - obs_limits[1, 3]

    obs_x_scale <- max(abs(obs_limits[, 1]))
    obs_y_scale <- max(abs(obs_limits[, 2]))
    obs_z_scale <- max(abs(obs_limits[, 3]))

    full_volume <- (obs_x_range * (1 + (2 * limit_scaler))) *
      (obs_y_range * (1 + (2 * limit_scaler))) *
      (obs_z_range * (1 + (2 * limit_scaler)))

    temp_data <- sim_mat[sim_mat[, 1] == time_points[time_idx], -1]
    data_limits <- colMinsMaxs(temp_data)

    x_range <- data_limits[2, 1] - data_limits[1, 1]
    y_range <- data_limits[2, 2] - data_limits[1, 2]
    z_range <- data_limits[2, 3] - data_limits[1, 3]

    x_min <- data_limits[1, 1] - (limit_scaler * x_range)
    x_max <- data_limits[2, 1] + (limit_scaler * x_range)
    y_min <- data_limits[1, 2] - (limit_scaler * y_range)
    y_max <- data_limits[2, 2] + (limit_scaler * y_range)
    z_min <- data_limits[1, 3] - (limit_scaler * z_range)
    z_max <- data_limits[2, 3] + (limit_scaler * z_range)

    x_candidates <- runif(n = n_points, min = x_min, max = x_max)
    y_candidates <- runif(n = n_points, min = y_min, max = y_max)
    z_candidates <- runif(n = n_points, min = z_min, max = z_max)

    candidates <- cbind(x_candidates, y_candidates, z_candidates)

    # Use partitioned models with interior identification
    temp_pme_part1 <- sim$models$pme_part$pme[[1]][[time_idx]]
    temp_pme_part2 <- sim$models$pme_part$pme[[2]][[time_idx]]
    temp_pme_embedding_part1 <- temp_pme_part1$embedding_map
    temp_pme_embedding_part2 <- temp_pme_part2$embedding_map
    temp_pme_coefs_part1 <- temp_pme_part1$coefs[[which.min(temp_pme_part1$MSD)]]
    temp_pme_coefs_part2 <- temp_pme_part2$coefs[[which.min(temp_pme_part2$MSD)]]
    temp_pme_params_part1 <- temp_pme_part1$parameterization[[which.min(temp_pme_part1$MSD)]]
    temp_pme_params_part2 <- temp_pme_part2$parameterization[[which.min(temp_pme_part2$MSD)]]

    lpme_est_part1 <- sim$models$lpme_part$lpme[[1]]
    lpme_est_part2 <- sim$models$lpme_part$lpme[[2]]

    temp_lpme_coefs_part1 <- lpme_est_part1$sol_coef_functions[[
      which.min(lpme_est_part1$msd)
    ]](time_points[time_idx])
    temp_lpme_coefs_part2 <- lpme_est_part2$sol_coef_functions[[
      which.min(lpme_est_part2$msd)
    ]](time_points[time_idx])

    temp_lpme_params_part1 <- lpme_est_part1$parameterization_list[[which.min(lpme_est_part1$msd)]]
    temp_lpme_params_part2 <- lpme_est_part2$parameterization_list[[which.min(lpme_est_part2$msd)]]
    temp_lpme_params_part1 <- temp_lpme_params_part1[
      temp_lpme_params_part1[, 1] == time_points[time_idx],
      -1
    ]
    temp_lpme_params_part2 <- temp_lpme_params_part2[
      temp_lpme_params_part2[, 1] == time_points[time_idx],
      -1
    ]

    lpme_n_knots_part1 <- nrow(temp_lpme_params_part1)
    lpme_n_knots_part2 <- nrow(temp_lpme_params_part2)
    d <- ncol(temp_lpme_params_part1)
    coef_mat_part1 <- matrix(
      temp_lpme_coefs_part1,
      lpme_n_knots_part1 + d + 1,
      byrow = TRUE
    )
    coef_mat_part2 <- matrix(
      temp_lpme_coefs_part2,
      lpme_n_knots_part2 + d + 1,
      byrow = TRUE
    )

    temp_lpme_embedding_part1 <- function(r) {
      t(coef_mat_part1[1:lpme_n_knots_part1, ]) %*% pme::etaFunc(r, temp_lpme_params_part1, 4 - d) +
        t(coef_mat_part1[(lpme_n_knots_part1 + 1):(lpme_n_knots_part1 + d + 1), ]) %*%
          matrix(c(1, r), ncol = 1)
    }

    temp_lpme_embedding_part2 <- function(r) {
      t(coef_mat_part2[1:lpme_n_knots_part2, ]) %*% pme::etaFunc(r, temp_lpme_params_part2, 4 - d) +
        t(coef_mat_part2[(lpme_n_knots_part2 + 1):(lpme_n_knots_part2 + d + 1), ]) %*%
          matrix(c(1, r), ncol = 1)
    }

    lpme_interior_part1 <- interior_identification(
      temp_lpme_embedding_part1,
      coef_mat_part1,
      temp_lpme_params_part1,
      candidates[candidates[, 1] > 0, ],
      c(0, 0, 0)
    )
    lpme_interior_part2 <- interior_identification(
      temp_lpme_embedding_part2,
      coef_mat_part2,
      temp_lpme_params_part2,
      candidates[candidates[, 1] <= 0, ],
      c(0, 0, 0)
    )

    lpme_interior <- c(lpme_interior_part1, lpme_interior_part2)

    lpme_interior_candidates <- rbind(
      candidates[candidates[, 1] > 0, ][lpme_interior_part1, ],
      candidates[candidates[, 1] <= 0, ][lpme_interior_part2, ]
    )
    
    pme_interior_part1 <- interior_identification(
      temp_pme_embedding_part1,
      temp_pme_coefs_part1,
      temp_pme_params_part1,
      candidates[candidates[, 1] > 0, ],
      c(0, 0, 0)
    )
    pme_interior_part2 <- interior_identification(
      temp_pme_embedding_part2,
      temp_pme_coefs_part2,
      temp_pme_params_part2,
      candidates[candidates[, 1] <= 0, ],
      c(0, 0, 0)
    )

    pme_interior <- c(pme_interior_part1, pme_interior_part2)
    pme_interior_candidates <- rbind(
      candidates[candidates[, 1] > 0, ][pme_interior_part1, ],
      candidates[candidates[, 1] <= 0, ][pme_interior_part2, ]
    )

    lpme_part_volumes[time_idx] <- mean(lpme_interior) * full_volume
    pme_part_volumes[time_idx] <- mean(pme_interior) * full_volume


    temp_pme <- sim$models$pme$pme[[time_idx]]
    temp_pme_params <- temp_pme$parameterization[[which.min(temp_pme$MSD)]]

    temp_lpme_params <- lpme_params[lpme_params[, 1] == time_points[time_idx], -1]

    pme_param_limits <- colMinsMaxs(temp_pme_params)
    lpme_param_limits <- colMinsMaxs(temp_lpme_params)

    pme_param_candidates1 <- runif(
      n = n_points, 
      min = pme_param_limits[1, 1], 
      max = pme_param_limits[2, 1]
    )
    pme_param_candidates2 <- runif(
      n = n_points, 
      min = pme_param_limits[1, 2], 
      max = pme_param_limits[2, 2]
    )
    lpme_param_candidates1 <- runif(
      n = n_points,
      min = lpme_param_limits[1, 1],
      max = lpme_param_limits[2, 1]
    )
    lpme_param_candidates2 <- runif(
      n = n_points,
      min = lpme_param_limits[1, 2],
      max = lpme_param_limits[2, 2]
    )


    temp_lpme_reconstructions <- map(
      seq_len(n_points),
      ~ lpme_est$embedding_map(
        c(
          time_points[time_idx], 
          lpme_param_candidates1[.x], 
          lpme_param_candidates2[.x]
        )
      )
    ) |>
      reduce(rbind)
    temp_lpme_reconstructions <- temp_lpme_reconstructions[, 2:4]

    temp_pme_reconstructions <- map(
      seq_len(n_points),
      ~ temp_pme$embedding_map(
        c(
          pme_param_candidates1[.x],
          pme_param_candidates2[.x]
        )
      )
    ) |>
      reduce(rbind)
    temp_pme_reconstructions <- temp_pme_reconstructions[, 1:3]

    temp_lpme_reconstructions[, 1] <- temp_lpme_reconstructions[, 1] * x_scale
    temp_lpme_reconstructions[, 2] <- temp_lpme_reconstructions[, 2] * y_scale
    temp_lpme_reconstructions[, 3] <- temp_lpme_reconstructions[, 3] * z_scale

    temp_pme_reconstructions[, 1] <- temp_pme_reconstructions[, 1] * x_scale
    temp_pme_reconstructions[, 2] <- temp_pme_reconstructions[, 2] * y_scale
    temp_pme_reconstructions[, 3] <- temp_pme_reconstructions[, 3] * z_scale

    temp_lpme_volume <- estimate_volume(
      round(temp_lpme_reconstructions, 1),
      voxel_volume = 0.1^3
    )
    temp_pme_volume <- estimate_volume(
      round(temp_pme_reconstructions, 1),
      voxel_volume = 0.1^3
    )

    lpme_volumes[time_idx] <- temp_lpme_volume$volume
    pme_volumes[time_idx] <- temp_pme_volume$volume
  }

  volume_list <- list(
    times = observed_times,
    observed = observed_volumes,
    true = true_volumes,
    lpme_augmented = lpme_volumes,
    pme_augmented = pme_volumes,
    lpme_part = lpme_part_volumes,
    pme_part = pme_part_volumes
  )

  return(volume_list)
}
