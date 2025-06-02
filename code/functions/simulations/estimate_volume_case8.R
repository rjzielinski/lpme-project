estimate_volume_case8 <- function(
  sim,
  time_trend_val,
  time_change_val,
  n_points,
  threshold = 0.005
) {
  require(alphashape3d, warn.conflicts = FALSE, quietly = TRUE)
  require(dplyr, warn.conflicts = FALSE, quietly = TRUE)
  require(geometry, warn.conflicts = FALSE, quietly = TRUE)
  require(Morpho, warn.conflicts = FALSE, quietly = TRUE)
  require(purrr, warn.conflicts = FALSE, quietly = TRUE)
  require(Rfast, warn.conflicts = FALSE, quietly = TRUE)
  require(Rvcg, warn.conflicts = FALSE, quietly = TRUE)

  source("functions/interior_identification.R")
  source("functions/estimate_volume_interior.R")
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

  sim_mat_pt1 <- sim_mat[sim_mat[, 2] > 0, ]
  sim_mat_pt2 <- sim_mat[sim_mat[, 2] <= 0, ]

  # sim_data_centers <- sim$processed_data$df_observed |>
  sim_data_centers <- sim$data$df |>
    group_by(time) |>
    summarize(
      max_x1 = max(abs(X1)),
      max_x2 = max(abs(X2)),
      max_x3 = max(abs(X3))
    )

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
    time_adjustments_scaled <- (time_adjustments / max(time_adjustments)) *
      time_change_val
  }
  # subtract from initial amplitude to replicate effects of atrophy
  true_radii <- 1 - time_adjustments_scaled
  true_volumes <- (4 / 3) * pi * true_radii^3

  ## Calcalate volumes from partitioned models

  lpme_volumes <- vector(mode = "numeric", length = length(time_points))
  pme_volumes <- vector(mode = "numeric", length = length(time_points))
  pc_volumes <- vector(mode = "numeric", length = length(time_points))

  lpme_part_volumes_mesh <- vector(
    mode = "numeric",
    length = length(time_points)
  )
  pme_part_volumes_mesh <- vector(
    mode = "numeric",
    length = length(time_points)
  )

  for (time_idx in seq_along(time_points)) {
    time_val <- time_points[time_idx]
    temp_lpme_reconstructions <- sim$models$lpme$reconstructions
    temp_lpme_reconstructions <- temp_lpme_reconstructions[
      temp_lpme_reconstructions[, 1] == time_val,
      2:4
    ]
    temp_lpme_reconstructions_scaled <- temp_lpme_reconstructions
    temp_lpme_reconstructions_scaled[, 1] <- temp_lpme_reconstructions_scaled[,
      1
    ] *
      as.numeric(sim_data_centers[time_idx, 2])
    temp_lpme_reconstructions_scaled[, 2] <- temp_lpme_reconstructions_scaled[,
      2
    ] *
      as.numeric(sim_data_centers[time_idx, 3])
    temp_lpme_reconstructions_scaled[, 3] <- temp_lpme_reconstructions_scaled[,
      3
    ] *
      as.numeric(sim_data_centers[time_idx, 4])

    temp_lpme_part_reconstructions <- sim$models$lpme_part$reconstructions
    temp_lpme_part_reconstructions <- temp_lpme_part_reconstructions[
      temp_lpme_part_reconstructions[, 1] == time_val,
      2:4
    ]
    temp_lpme_part_reconstructions_scaled <- temp_lpme_part_reconstructions
    temp_lpme_part_reconstructions_scaled[,
      1
    ] <- temp_lpme_part_reconstructions_scaled[,
      1
    ] *
      as.numeric(sim_data_centers[time_idx, 2])
    temp_lpme_part_reconstructions_scaled[,
      2
    ] <- temp_lpme_part_reconstructions_scaled[,
      2
    ] *
      as.numeric(sim_data_centers[time_idx, 3])
    temp_lpme_part_reconstructions_scaled[,
      3
    ] <- temp_lpme_part_reconstructions_scaled[,
      3
    ] *
      as.numeric(sim_data_centers[time_idx, 4])

    temp_pme_reconstructions <- sim$models$pme$reconstructions
    temp_pme_reconstructions <- temp_pme_reconstructions[
      temp_pme_reconstructions[, 1] == time_val,
      2:4
    ]
    temp_pme_reconstructions_scaled <- temp_pme_reconstructions
    temp_pme_reconstructions_scaled[, 1] <- temp_pme_reconstructions_scaled[,
      1
    ] *
      as.numeric(sim_data_centers[time_idx, 2])
    temp_pme_reconstructions_scaled[, 2] <- temp_pme_reconstructions_scaled[,
      2
    ] *
      as.numeric(sim_data_centers[time_idx, 3])
    temp_pme_reconstructions_scaled[, 3] <- temp_pme_reconstructions_scaled[,
      3
    ] *
      as.numeric(sim_data_centers[time_idx, 4])

    temp_pme_part_reconstructions <- sim$models$pme_part$reconstructions
    temp_pme_part_reconstructions <- temp_pme_part_reconstructions[
      temp_pme_part_reconstructions[, 1] == time_val,
      2:4
    ]
    temp_pme_part_reconstructions_scaled <- temp_pme_part_reconstructions
    temp_pme_part_reconstructions_scaled[,
      1
    ] <- temp_pme_part_reconstructions_scaled[,
      1
    ] *
      as.numeric(sim_data_centers[time_idx, 2])
    temp_pme_part_reconstructions_scaled[,
      2
    ] <- temp_pme_part_reconstructions_scaled[,
      2
    ] *
      as.numeric(sim_data_centers[time_idx, 3])
    temp_pme_part_reconstructions_scaled[,
      3
    ] <- temp_pme_part_reconstructions_scaled[,
      3
    ] *
      as.numeric(sim_data_centers[time_idx, 4])

    temp_pc_reconstructions <- sim$models$pc_part$reconstructions
    temp_pc_reconstructions <- temp_pc_reconstructions[
      temp_pc_reconstructions[, 1] == time_val,
      2:4
    ]
    temp_pc_reconstructions_scaled <- temp_pc_reconstructions
    temp_pc_reconstructions_scaled[, 1] <- temp_pc_reconstructions_scaled[,
      1
    ] *
      as.numeric(sim_data_centers[time_idx, 2])
    temp_pc_reconstructions_scaled[, 2] <- temp_pc_reconstructions_scaled[,
      2
    ] *
      as.numeric(sim_data_centers[time_idx, 3])
    temp_pc_reconstructions_scaled[, 3] <- temp_pc_reconstructions_scaled[,
      3
    ] *
      as.numeric(sim_data_centers[time_idx, 4])

    # test a wide range of alpha values
    temp_lpme_ashape <- ashape3d(
      temp_lpme_reconstructions_scaled,
      alpha = exp(seq(-3, 2, 0.1))
    )
    temp_lpme_volumes <- volume_ashape3d(
      temp_lpme_ashape,
      indexAlpha = "all"
    )

    # volumes tend to increase drastically as alpha first increases
    # then volumes plateau - find where they stabilize
    temp_lpme_volume_change <- c(
      NA,
      ((lead(temp_lpme_volumes) - temp_lpme_volumes) / temp_lpme_volumes)[
        -length(temp_lpme_volumes)
      ]
    )
    # once volume has stabilized, take the median volume
    lpme_volumes[time_idx] <- median(
      temp_lpme_volumes[temp_lpme_volume_change < threshold],
      na.rm = TRUE
    )

    temp_lpme_part_ashape <- ashape3d(
      temp_lpme_part_reconstructions_scaled,
      alpha = exp(seq(-3, 2, 0.1))
    )
    temp_lpme_part_volumes <- volume_ashape3d(
      temp_lpme_part_ashape,
      indexAlpha = "all"
    )
    temp_lpme_part_volume_change <- c(
      NA,
      ((lead(temp_lpme_part_volumes) - temp_lpme_part_volumes) /
        temp_lpme_part_volumes)[-length(temp_lpme_part_volumes)]
    )
    lpme_part_volumes_mesh[time_idx] <- median(
      temp_lpme_part_volumes[temp_lpme_part_volume_change < threshold],
      na.rm = TRUE
    )

    temp_pme_ashape <- ashape3d(
      temp_pme_reconstructions_scaled,
      alpha = exp(seq(-3, 2, 0.1))
    )
    temp_pme_volumes <- volume_ashape3d(
      temp_pme_ashape,
      indexAlpha = "all"
    )
    temp_pme_volume_change <- c(
      NA,
      ((lead(temp_pme_volumes) - temp_pme_volumes) / temp_pme_volumes)[
        -length(temp_pme_volumes)
      ]
    )
    pme_volumes[time_idx] <- median(
      temp_pme_volumes[temp_pme_volume_change < threshold],
      na.rm = TRUE
    )

    temp_pme_part_ashape <- ashape3d(
      temp_pme_part_reconstructions_scaled,
      alpha = exp(seq(-3, 2, 0.1))
    )
    temp_pme_part_volumes <- volume_ashape3d(
      temp_pme_part_ashape,
      indexAlpha = "all"
    )
    temp_pme_part_volume_change <- c(
      NA,
      ((lead(temp_pme_part_volumes) - temp_pme_part_volumes) /
        temp_pme_part_volumes)[-length(temp_pme_part_volumes)]
    )
    pme_part_volumes_mesh[time_idx] <- median(
      temp_pme_part_volumes[temp_pme_part_volume_change < threshold],
      na.rm = TRUE
    )

    temp_pc_ashape <- ashape3d(
      temp_pc_reconstructions_scaled,
      alpha = exp(seq(-3, 2, 0.1))
    )
    temp_pc_volumes <- volume_ashape3d(
      temp_pc_ashape,
      indexAlpha = "all"
    )
    temp_pc_volume_change <- c(
      NA,
      ((lead(temp_pc_volumes) - temp_pc_volumes) / temp_pc_volumes)[
        -length(temp_pc_volumes)
      ]
    )
    pc_volumes[time_idx] <- median(
      temp_pc_volumes[temp_pc_volume_change < threshold],
      na.rm = TRUE
    )
  }

  lpme_part_volumes <- estimate_volume_interior_lpme(
    sim$models$lpme_part$lpme,
    list(sim_mat_pt1, sim_mat_pt2),
    time_points,
    n_points = 10000,
    data_max = sim_data_centers,
    limit_scaler = 0.05,
    partition_index = 1
  )

  pme_part_volumes <- estimate_volume_interior_pme(
    sim$models$pme_part$pme,
    list(sim_mat_pt1, sim_mat_pt2),
    time_points,
    n_points = 10000,
    data_max = sim_data_centers,
    limit_scaler = 0.05,
    partition_index = 1
  )

  volume_list <- list(
    times = observed_times,
    observed = observed_volumes,
    true = true_volumes,
    lpme_augmented = lpme_volumes,
    pme_augmented = pme_volumes,
    lpme_part = lpme_part_volumes$volumes,
    pme_part = pme_part_volumes$volumes,
    lpme_part_mesh = lpme_part_volumes_mesh,
    pme_part_mesh = pme_part_volumes_mesh,
    pc_part = pc_volumes
  )

  return(volume_list)
}
