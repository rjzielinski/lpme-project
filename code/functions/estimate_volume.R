estimate_volume <- function(data, voxel_volume) {
  require(Rfast, quietly = TRUE, warn.conflicts = FALSE)

  indices <- matrix(ncol = 3)
  data_limits <- colMinsMaxs(data)

  x_seq <- seq(data_limits[1, 1], data_limits[2, 1], by = 0.1)
  y_seq <- seq(data_limits[1, 2], data_limits[2, 2], by = 0.1)
  z_seq <- seq(data_limits[1, 3], data_limits[2, 3], by = 0.1)

  fill_candidates <- expand.grid(x_seq, y_seq, z_seq)
  interior_point <- vector(mode = "logical", length = nrow(fill_candidates))

  for (row_idx in seq_len(nrow(fill_candidates))) {
    interior_point[row_idx] <- check_interior(
      unlist(fill_candidates[row_idx, ]),
      data
    )
  }

  interior_voxels <- fill_candidates[interior_point, ]
  est_volume <- sum(interior_point) * voxel_volume
  volume_out <- list(
    volume = est_volume,
    interior_voxels = interior_voxels
  )
  return(volume_out)
}

check_interior <- function(x, data) {
  data_subsets <- list(
    matrix(
      data[near(data[, 2], x[2]) & near(data[, 3], x[3]), ],
      ncol = ncol(data)
    ),
    matrix(
      data[near(data[, 1], x[1]) & near(data[, 3], x[3]), ],
      ncol = ncol(data)
    ),
    matrix(
      data[near(data[, 1], x[1]) & near(data[, 2], x[2]), ],
      ncol = ncol(data)
    )
  )

  subset_dims <- map(
    data_subsets,
    ~ nrow(.x)
  ) |>
    reduce(c)

  interior <- vector(mode = "logical", length = length(subset_dims))

  for (dim_idx in seq_along(data_subsets)) {
    if (subset_dims[dim_idx] > 0) {
      if (subset_dims[dim_idx] > 1) {
        min_val <- min(data_subsets[[dim_idx]][, dim_idx])
        max_val <- max(data_subsets[[dim_idx]][, dim_idx])
        if ((x[dim_idx] >= min_val) && (x[dim_idx] <= max_val)) {
          interior[dim_idx] <- TRUE
        }
      } else {
        val <- max(abs(data_subsets[[dim_idx]][, dim_idx]))
        if (abs(x[dim_idx]) <= abs(val)) {
          interior[dim_idx] <- TRUE
        }
      }
    } else {
      interior[dim_idx] <- FALSE
    }
  }

  if (sum(interior) >= 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


