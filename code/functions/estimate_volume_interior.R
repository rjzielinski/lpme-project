estimate_volume_interior_lpme <- function(
  lpme_list,
  data,
  time_points,
  n_points,
  data_max,
  limit_scaler,
  partition_index,
  augment = FALSE
) {
  volumes <- vector(mode = "numeric", length = length(time_points))

  voxels <- matrix(data = NA, ncol = 4)

  lpme_est_pt1 <- lpme_list[[1]]
  lpme_est_pt2 <- lpme_list[[2]]

  for (time_idx in seq_along(time_points)) {
    x_scale <- as.numeric(data_max[time_idx, 2])
    y_scale <- as.numeric(data_max[time_idx, 3])
    z_scale <- as.numeric(data_max[time_idx, 4])

    full_volume <- (2 * (x_scale * (1 + limit_scaler))) *
      (2 * (y_scale * (1 + limit_scaler))) *
      (2 * (z_scale * (1 + limit_scaler)))

    x_candidates <- runif(
      n = n_points,
      min = -x_scale * (1 + limit_scaler),
      max = x_scale * (1 + limit_scaler)
    )

    y_candidates <- runif(
      n = n_points,
      min = -y_scale * (1 + limit_scaler),
      max = y_scale * (1 + limit_scaler)
    )

    z_candidates <- runif(
      n = n_points,
      min = -z_scale * (1 + limit_scaler),
      max = z_scale * (1 + limit_scaler)
    )

    candidates <- cbind(x_candidates, y_candidates, z_candidates)
    candidate_means <- c(
      mean(candidates[, 1]),
      mean(candidates[, 2]),
      mean(candidates[, 3])
    )
    candidate_maxs <- c(
      max(abs(candidates[, 1] - candidate_means[1])),
      max(abs(candidates[, 2] - candidate_means[2])),
      max(abs(candidates[, 3] - candidate_means[3]))
    )

    candidates_scaled <- candidates
    candidates_scaled[, 1] <- (candidates_scaled[, 1] - candidate_means[1]) /
      candidate_maxs[1]
    candidates_scaled[, 2] <- (candidates_scaled[, 2] - candidate_means[2]) /
      candidate_maxs[2]
    candidates_scaled[, 3] <- (candidates_scaled[, 3] - candidate_means[3]) /
      candidate_maxs[3]

    temp_lpme_coefs_pt1 <- lpme_est_pt1$sol_coef_functions[[which.min(
      lpme_est_pt1$msd
    )]](time_points[time_idx])
    temp_lpme_coefs_pt2 <- lpme_est_pt2$sol_coef_functions[[which.min(
      lpme_est_pt2$msd
    )]](time_points[time_idx])
    temp_lpme_params_pt1 <- lpme_est_pt1$parameterization_list[[which.min(
      lpme_est_pt1$msd
    )]]
    temp_lpme_params_pt1 <- temp_lpme_params_pt1[
      temp_lpme_params_pt1[, 1] == time_points[time_idx],
      -1
    ]
    temp_lpme_params_pt2 <- lpme_est_pt2$parameterization_list[[which.min(
      lpme_est_pt2$msd
    )]]
    temp_lpme_params_pt2 <- temp_lpme_params_pt2[
      temp_lpme_params_pt2[, 1] == time_points[time_idx],
      -1
    ]

    lpme_n_knots_pt1 <- nrow(temp_lpme_params_pt1)
    lpme_n_knots_pt2 <- nrow(temp_lpme_params_pt2)
    d <- ncol(temp_lpme_params_pt1)

    coef_mat_pt1 <- matrix(
      temp_lpme_coefs_pt1,
      lpme_n_knots_pt1 + d + 1,
      byrow = TRUE
    )
    coef_mat_pt2 <- matrix(
      temp_lpme_coefs_pt2,
      lpme_n_knots_pt2 + d + 1,
      byrow = TRUE
    )

    temp_lpme_embedding_pt1 <- function(r) {
      t(coef_mat_pt1[1:lpme_n_knots_pt1, ]) %*%
        pme::etaFunc(r, temp_lpme_params_pt1, 4 - d) +
        t(coef_mat_pt1[(lpme_n_knots_pt1 + 1):(lpme_n_knots_pt1 + d + 1), ]) %*%
          matrix(c(1, r), ncol = 1)
    }

    temp_lpme_embedding_pt2 <- function(r) {
      t(coef_mat_pt2[1:lpme_n_knots_pt2, ]) %*%
        pme::etaFunc(r, temp_lpme_params_pt2, 4 - d) +
        t(coef_mat_pt2[(lpme_n_knots_pt2 + 1):(lpme_n_knots_pt2 + d + 1), ]) %*%
          matrix(c(1, r), ncol = 1)
    }

    lpme_interior_pt1 <- interior_identification(
      temp_lpme_embedding_pt1,
      coef_mat_pt1,
      temp_lpme_params_pt1,
      candidates_scaled[candidates_scaled[, partition_index] > 0, ],
      c(0, 0, 0)
    )
    lpme_interior_pt2 <- interior_identification(
      temp_lpme_embedding_pt2,
      coef_mat_pt2,
      temp_lpme_params_pt2,
      candidates_scaled[candidates_scaled[, partition_index] <= 0, ],
      c(0, 0, 0)
    )

    lpme_interior <- c(lpme_interior_pt1, lpme_interior_pt2)
    lpme_interior_candidates <- rbind(
      candidates[candidates_scaled[, partition_index] > 0, ][
        lpme_interior_pt1,
      ],
      candidates[candidates_scaled[, partition_index] <= 0, ][
        lpme_interior_pt2,
      ]
    )

    voxels <- rbind(
      voxels,
      cbind(time_points[time_idx], lpme_interior_candidates)
    )
    volumes[time_idx] <- mean(lpme_interior) * full_volume
  }

  volume_out <- list(
    volumes = volumes,
    interior_voxels = voxels
  )
}

estimate_volume_interior_pme <- function(
  pme_list,
  data,
  time_points,
  n_points,
  data_max,
  limit_scaler,
  augment = FALSE,
  partition_index
) {
  volumes <- vector(mode = "numeric", length = length(time_points))

  voxels <- matrix(data = NA, ncol = 4)

  pme_est_pt1 <- pme_list[[1]]
  pme_est_pt2 <- pme_list[[2]]

  for (time_idx in seq_along(time_points)) {
    x_scale <- as.numeric(data_max[time_idx, 2])
    y_scale <- as.numeric(data_max[time_idx, 3])
    z_scale <- as.numeric(data_max[time_idx, 4])

    full_volume <- (2 * (x_scale * (1 + limit_scaler))) *
      (2 * (y_scale * (1 + limit_scaler))) *
      (2 * (z_scale * (1 + limit_scaler)))

    x_candidates <- runif(
      n = n_points,
      min = -x_scale * (1 + limit_scaler),
      max = x_scale * (1 + limit_scaler)
    )

    y_candidates <- runif(
      n = n_points,
      min = -y_scale * (1 + limit_scaler),
      max = y_scale * (1 + limit_scaler)
    )

    z_candidates <- runif(
      n = n_points,
      min = -z_scale * (1 + limit_scaler),
      max = z_scale * (1 + limit_scaler)
    )

    candidates <- cbind(x_candidates, y_candidates, z_candidates)
    candidate_means <- c(
      mean(candidates[, 1]),
      mean(candidates[, 2]),
      mean(candidates[, 3])
    )
    candidate_maxs <- c(
      max(abs(candidates[, 1] - candidate_means[1])),
      max(abs(candidates[, 2] - candidate_means[2])),
      max(abs(candidates[, 3] - candidate_means[3]))
    )

    candidates_scaled <- candidates
    candidates_scaled[, 1] <- (candidates_scaled[, 1] - candidate_means[1]) /
      candidate_maxs[1]
    candidates_scaled[, 2] <- (candidates_scaled[, 2] - candidate_means[2]) /
      candidate_maxs[2]
    candidates_scaled[, 3] <- (candidates_scaled[, 3] - candidate_means[3]) /
      candidate_maxs[3]

    temp_pme_pt1 <- pme_est_pt1[[time_idx]]
    temp_pme_pt2 <- pme_est_pt2[[time_idx]]
    temp_pme_embedding_pt1 <- temp_pme_pt1$embedding_map
    temp_pme_embedding_pt2 <- temp_pme_pt2$embedding_map
    temp_pme_coefs_pt1 <- temp_pme_pt1$coefs[[which.min(temp_pme_pt1$MSD)]]
    temp_pme_coefs_pt2 <- temp_pme_pt2$coefs[[which.min(temp_pme_pt2$MSD)]]
    temp_pme_params_pt1 <- temp_pme_pt1$parameterization[[which.min(
      temp_pme_pt1$MSD
    )]]
    temp_pme_params_pt2 <- temp_pme_pt2$parameterization[[which.min(
      temp_pme_pt2$MSD
    )]]

    pme_interior_pt1 <- interior_identification(
      temp_pme_embedding_pt1,
      temp_pme_coefs_pt1,
      temp_pme_params_pt1,
      candidates_scaled[candidates_scaled[, partition_index] > 0, ],
      c(0, 0, 0)
    )
    pme_interior_pt2 <- interior_identification(
      temp_pme_embedding_pt2,
      temp_pme_coefs_pt2,
      temp_pme_params_pt2,
      candidates_scaled[candidates_scaled[, partition_index] <= 0, ],
      c(0, 0, 0)
    )

    pme_interior <- c(pme_interior_pt1, pme_interior_pt2)
    pme_interior_candidates <- rbind(
      candidates[candidates_scaled[, partition_index] > 0, ][
        pme_interior_pt1,
      ],
      candidates[candidates_scaled[, partition_index] <= 0, ][
        pme_interior_pt2,
      ]
    )

    voxels <- rbind(
      voxels,
      cbind(time_points[time_idx], pme_interior_candidates)
    )
    volumes[time_idx] <- mean(pme_interior) * full_volume
  }

  volume_out <- list(
    volumes = volumes,
    interior_voxels = voxels
  )
}
