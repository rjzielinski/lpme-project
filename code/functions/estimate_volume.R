estimate_volume <- function(data, voxel_volume) {
  require(Rfast, quietly = TRUE, warn.conflicts = FALSE)

  indices <- matrix(ncol = 3)
  data_limits <- colMinsMaxs(data)

  x_seq <- seq(data_limits[1, 1], data_limits[2, 1])
  y_seq <- seq(data_limits[1, 2], data_limits[2, 2])
  z_seq <- seq(data_limits[1, 3], data_limits[2, 3])

  # fill_candidates <- expand.grid(x_seq, y_seq, z_seq)
  fill_candidates <- matrix(ncol = 3)
  x_mins <- matrix(nrow = length(y_seq), ncol = length(z_seq))
  x_maxs <- matrix(nrow = length(y_seq), ncol = length(z_seq))

  for (y_idx in seq_along(y_seq)) {
    for (z_idx in seq_along(z_seq)) {
      y_val <- y_seq[y_idx]
      z_val <- z_seq[z_idx]
      data_red <- data[data[, 2] == y_val & data[, 3] == z_val, ] |>
        matrix(ncol = ncol(data))

      if (nrow(data_red) == 0) {
        x_mins[y_idx, z_idx] <- NA
        x_maxs[y_idx, z_idx] <- NA
      } else if (nrow(data_red) > 1) {
        min_val <- min(unique(data_red[, 1]))
        max_val <- max(unique(data_red[, 1]))
        if (max_val - min_val > 1) {
          x_mins[y_idx, z_idx] <- min(unique(data_red[, 1]))
          x_maxs[y_idx, z_idx] <- max(unique(data_red[, 1]))
        }
      } else {
        x_mins[y_idx, z_idx] <- NA
        x_maxs[y_idx, z_idx] <- NA
      }
    }
  }

  for (y_idx in seq_along(y_seq)) {
    for (z_idx in seq_along(z_seq)) {
      y_val <- y_seq[y_idx]
      z_val <- z_seq[z_idx]
      data_red <- data[data[, 2] == y_val & data[, 3] == z_val, ] |>
        matrix(ncol = ncol(data))

      if (nrow(data_red) == 0) {
        next
      } else if (nrow(data_red) > 1) {
        next
      } else {
        x_min_vals <- get_surrounding_vals(x_mins, y_idx, z_idx)
        x_max_vals <- get_surrounding_vals(x_maxs, y_idx, z_idx)
        if (sum(is.na(x_min_vals)) != length(x_min_vals)) {
          total_min <- min(x_min_vals, na.rm = TRUE)
        } else {
          total_min <- NA
          next
        }
        if (sum(is.na(x_max_vals)) != length(x_max_vals)) {
          total_max <- max(x_max_vals, na.rm = TRUE)
        } else {
          total_max <- NA
          next
        }
        # min_mean <- mean(x_min_vals, na.rm = TRUE)
        # max_mean <- mean(x_max_vals, na.rm = TRUE)
        if (abs(data_red[, 1] - total_min) > abs(data_red[, 1] - total_max)) {
          x_maxs[y_idx, z_idx] <- data_red[, 1]
        } else {
          x_mins[y_idx, z_idx] <- data_red[, 1]
        }
      }
    }
  }

  missing_present <- TRUE
  n_missing_prev <- nrow(x_mins) * ncol(x_mins)
  missingness_thresholds <- c(0.25, 0.33, 0.5, 0.66, 0.75, 0.99)
  threshold_idx <- 1

  while (missing_present == TRUE) {
    missing_mins <- matrix(nrow = nrow(x_mins), ncol = ncol(x_mins))
    missing_maxs <- matrix(nrow = nrow(x_maxs), ncol = ncol(x_maxs))
    n_missing <- 0
    for (y_idx in seq_along(y_seq)) {
      for (z_idx in seq_along(z_seq)) {
        if (is.na(x_mins[y_idx, z_idx])) {
          data_red <- data[
            data[, 2] == y_seq[y_idx] & data[, 3] == z_seq[z_idx],
          ] |>
            matrix(ncol = ncol(data))
          if (nrow(data_red) > 0) {
            n_missing <- n_missing + 1
            missing_mins[y_idx, z_idx] <- TRUE
          } else if (!is.na(x_maxs[y_idx, z_idx])) {
            n_missing <- n_missing + 1
            missing_mins[y_idx, z_idx] <- TRUE
          } else {
            n_present_y_lower <- x_mins[y_idx, 1:z_idx]
            n_present_y_upper <- x_mins[y_idx, z_idx:length(z_seq)]
            if ((sum(!is.na(n_present_y_lower)) > 0) && (sum(!is.na(n_present_y_upper)) > 0)) {
              n_missing <- n_missing + 1
              missing_mins[y_idx, z_idx] <- TRUE
            } else {
              n_present_z_lower <- x_mins[1:y_idx, z_idx]
              n_present_z_upper <- x_mins[y_idx:length(y_seq), z_idx]
              if ((sum(!is.na(n_present_z_lower)) > 0) && (sum(!is.na(n_present_z_upper)) > 0)) {
                n_missing <- n_missing + 1
                missing_mins[y_idx, z_idx] <- TRUE
              } else {
                missing_mins[y_idx, z_idx] <- FALSE
              }
            }
          }
        } else {
          missing_mins[y_idx, z_idx] <- FALSE
        }
        #   x_min_vals <- get_surrounding_vals(x_mins, y_idx, z_idx)
        #   num_missing <- sum(is.na(x_min_vals))
        #   if ((num_missing < 4) | (length(x_min_vals) - num_missing) >= 5) {
        #     n_missing <- n_missing + 1
        #     missing_mins[y_idx, z_idx] <- TRUE
        #   } else {
        #     missing_mins[y_idx, z_idx] <- FALSE
        #   }
        # } else {
        #   missing_mins[y_idx, z_idx] <- FALSE
        # }
        # if (is.na(x_maxs[y_idx, z_idx])) {
        #   x_max_vals <- get_surrounding_vals(x_maxs, y_idx, z_idx)
        #   num_missing <- sum(is.na(x_max_vals))
        #   if ((num_missing < 4) | (length(x_max_vals) - num_missing) >= 5) {
        #     n_missing <- n_missing + 1
        #     missing_maxs[y_idx, z_idx] <- TRUE
        #   } else {
        #     missing_maxs[y_idx, z_idx] <- FALSE
        #   }
        # } else {
        #   missing_maxs[y_idx, z_idx] <- FALSE
        # }
        if (is.na(x_maxs[y_idx, z_idx])) {
          data_red <- data[data[, 2] == y_seq[y_idx] & data[, 3] == z_seq[z_idx], ] |>
            matrix(ncol = ncol(data))
          if (nrow(data_red) > 0) {
            n_missing <- n_missing + 1
            missing_maxs[y_idx, z_idx] <- TRUE
          } else if (!is.na(x_mins[y_idx, z_idx])) {
            n_missing <- n_missing + 1
            missing_maxs[y_idx, z_idx] <- TRUE
          } else {
            n_present_y_lower <- x_maxs[y_idx, 1:z_idx]
            n_present_y_upper <- x_maxs[y_idx, z_idx:length(z_seq)]
            if ((sum(!is.na(n_present_y_lower)) > 0) && (sum(!is.na(n_present_y_upper)) > 0)) {
              n_missing <- n_missing + 1
              missing_maxs[y_idx, z_idx] <- TRUE
            } else {
              n_present_z_lower <- x_maxs[1:y_idx, z_idx]
              n_present_z_upper <- x_maxs[y_idx:length(y_seq), z_idx]
              if ((sum(!is.na(n_present_z_lower)) > 0) && (sum(!is.na(n_present_z_upper)) > 0)) {
                n_missing <- n_missing + 1
                missing_maxs[y_idx, z_idx] <- TRUE
              } else {
                missing_maxs[y_idx, z_idx] <- FALSE
              }
            }
          }
        } else {
          missing_maxs[y_idx, z_idx] <- FALSE
        }
      }
    }

    if (n_missing == 0) {
      missing_present <- FALSE
    } else {
      missing_present <- TRUE
    }

    if (!(n_missing < n_missing_prev)) {
      if (threshold_idx < length(missingness_thresholds)) {
        threshold_idx <- threshold_idx + 1
      } else {
        break
      }
    }
    n_missing_prev <- n_missing

    for (y_idx in seq_along(y_seq)) {
      for (z_idx in seq_along(z_seq)) {
        if (missing_mins[y_idx, z_idx] == TRUE) {
          x_min_vals <- get_surrounding_vals(x_mins, y_idx, z_idx)
          if (sum(is.na(x_min_vals)) / length(x_min_vals) <= missingness_thresholds[threshold_idx]) {
            # if (sum(!is.na(x_min_vals)) > 1) {
            # x_mins[y_idx, z_idx] <- min(x_min_vals, na.rm = TRUE)
            x_mins[y_idx, z_idx] <- round(mean(x_min_vals, na.rm = TRUE))
            missing_mins[y_idx, z_idx] <- FALSE
            # n_missing <- n_missing - 1
          }
        }
        if (missing_maxs[y_idx, z_idx] == TRUE) {
          x_max_vals <- get_surrounding_vals(x_maxs, y_idx, z_idx)
          if (sum(!is.na(x_max_vals)) / length(x_max_vals) <= missingness_thresholds[threshold_idx]) {
            # if (sum(!is.na(x_max_vals)) > 1) {
            # x_maxs[y_idx, z_idx] <- max(x_max_vals, na.rm = TRUE)
            x_maxs[y_idx, z_idx] <- round(mean(x_max_vals, na.rm = TRUE))
            missing_maxs[y_idx, z_idx] <- FALSE
            # n_missing <- n_missing - 1
          }
        }
      }
    }

    # n_missing <- 0
    # for (y_idx in 1:length(y_seq)) {
    #   for (z_idx in 1:length(z_seq)) {
    #     if (is.na(x_mins[y_idx, z_idx])) {
    #       x_min_vals <- get_surrounding_vals(x_mins, y_idx, z_idx)
    #       num_missing <- sum(is.na(x_min_vals))
    #       if ((num_missing < 4) | (length(x_min_vals) - num_missing) >= 5) {
    #         n_missing <- n_missing + 1
    #         missing_mins[y_idx, z_idx] <- TRUE
    #       } else {
    #         missing_mins[y_idx, z_idx] <- FALSE
    #       }
    #     } else {
    #       missing_mins[y_idx, z_idx] <- FALSE
    #     }
    #     if (is.na(x_maxs[y_idx, z_idx])) {
    #       x_max_vals <- get_surrounding_vals(x_maxs, y_idx, z_idx)
    #       num_missing <- sum(is.na(x_max_vals))
    #       if ((num_missing < 4) | (length(x_max_vals) - num_missing) >= 5) {
    #         n_missing <- n_missing + 1
    #         missing_maxs[y_idx, z_idx] <- TRUE
    #       } else {
    #         missing_maxs[y_idx, z_idx] <- FALSE
    #       }
    #     } else {
    #       missing_maxs[y_idx, z_idx] <- FALSE
    #     }
    #   }
    # }

    # if (n_missing > 0) {
    #   missing_present <- TRUE
    # } else {
    #   missing_present <- FALSE
    # }
    # print(sum(!is.na(x_mins)))
    # print(sum(!is.na(x_maxs)))
  }

  volume <- 0
  for (y_idx in seq_along(y_seq)) {
    for (z_idx in seq_along(z_seq)) {
      x_min_val <- x_mins[y_idx, z_idx]
      x_max_val <- x_maxs[y_idx, z_idx]

      if (!is.na(x_min_val) && !is.na((x_max_val))) {
        x_vals <- seq(min(x_min_val, x_max_val), max(x_min_val, x_max_val), by = 1)
        for (x_val in x_vals) {
          fill_candidates <- rbind(
            fill_candidates,
            c(x_val, y_seq[y_idx], z_seq[z_idx])
          )
        }
      }
    }
  }
  fill_candidates <- fill_candidates[-1, ]
  volume <- nrow(fill_candidates) * voxel_volume
  return(list(volume, fill_candidates))
}

get_surrounding_vals <- function(mat, x_idx, y_idx) {
  if (x_idx == 1) {
    if (y_idx == 1) {
      out <- mat[x_idx:(x_idx + 1), y_idx:(y_idx + 1)]
    } else if (y_idx == ncol(mat)) {
      out <- mat[x_idx:(x_idx + 1), (y_idx - 1):y_idx]
    } else {
      out <- mat[x_idx:(x_idx + 1), (y_idx - 1):(y_idx + 1)]
    }
  } else if (x_idx == nrow(mat)) {
    if (y_idx == 1) {
      out <- mat[(x_idx - 1):x_idx, y_idx:(y_idx + 1)]
    } else if (y_idx == ncol(mat)) {
      out <- mat[(x_idx - 1):x_idx, (y_idx - 1):y_idx]
    } else {
      out <- mat[(x_idx - 1):x_idx, (y_idx - 1):(y_idx + 1)]
    }
  } else {
    if (y_idx == 1) {
      out <- mat[(x_idx - 1):(x_idx + 1), y_idx:(y_idx + 1)]
    } else if (y_idx == ncol(mat)) {
      out <- mat[(x_idx - 1):(x_idx + 1), (y_idx - 1):y_idx]
    } else {
      out <- mat[(x_idx - 1):(x_idx + 1), (y_idx - 1):(y_idx + 1)]
    }
  }
  out
}
