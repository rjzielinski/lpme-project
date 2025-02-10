fit_models <- function(data, case, d, D) {
  ##### MANIFOLD LEARNING #####

  require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  require(pme, quietly = TRUE, warn.conflicts = FALSE)
  require(princurve, quietly = TRUE, warn.conflicts = FALSE)
  require(purrr, quietly = TRUE, warn.conflicts = FALSE)

  if (case != 8) {
    mat <- as.matrix(data$df_observed)

    lpme_results <- list()

    lpme_results[[1]] <- lpme(
      mat,
      d,
      verbose = FALSE,
      print_plots = FALSE
    )
    lpme_reconstructions <- calculate_lpme_reconstructions(
      lpme_results[[1]],
      mat
    )

    # fit PME and appropriate principal curve/surface algorithms and
    # calculate reconstructions for each
    pme_result_list <- list()
    pc_result_list <- list()

    pme_reconstruction_list <- list()
    pc_reconstruction_list <- list()

    smoothing_options <- c("smooth_spline", "lowess", "periodic_lowess")

    time_values <- unique(data$df$time)
    for (time_idx in seq_along(time_values)) {
      temp_data <- mat[mat[, 1] == time_values[time_idx], -1]
      pme_result_list[[time_idx]] <- pme(temp_data, d = d)
      pme_reconstruction_list[[time_idx]] <- calculate_pme_reconstructions(
        pme_result_list[[time_idx]],
        temp_data
      )
      pme_reconstruction_list[[time_idx]] <- cbind(
        time_values[time_idx],
        pme_reconstruction_list[[time_idx]]
      )
      if (d == 1) {
        principal_curves <- list()
        pc_error <- vector()
        for (smoother_idx in seq_along(smoothing_options)) {
          principal_curves[[smoother_idx]] <- principal_curve(
            temp_data,
            smoother = smoothing_options[[smoother_idx]]
          )
          pc_error[smoother_idx] <- principal_curves[[smoother_idx]]$dist
        }
        opt_principal_curve <- which.min(pc_error)
        pc_result_list[[time_idx]] <- principal_curves[[opt_principal_curve]]
        pc_reconstruction_list[[time_idx]] <- cbind(
          time_values[time_idx],
          pc_result_list[[time_idx]]$s
        )
      } else if (d == 2) {
        if (dim(temp_data)[2] == 3) {
          principal_surface <- prinSurf(temp_data)
          surface_mse <- map(
            seq_along(principal_surface),
            ~ principal_surface[[.x]]$MSE
          ) |>
            unlist()
          opt_surface <- which.min(surface_mse)
          # using opt_surface + 2 as index because first two list entries are NULL
          pc_result_list[[time_idx]] <- principal_surface[[opt_surface + 2]]
          pc_reconstruction_list[[time_idx]] <- cbind(
            time_values[time_idx],
            principal_surface[[opt_surface + 2]]$PS
          )
        } else if (dim(temp_data)[2] > 3) {
          # As coded, the principal surface function assumes that D = 3
          # This is not the case when considering augmented data
          pc_result_list[[time_idx]] <- NULL
          pc_reconstruction_list[[time_idx]] <- NULL
        }
      }
    }
  }

  if (case == 8) {
    mat_part1 <- data$df_observed |>
      filter(X1 > 0) |>
      select(time, contains("X")) |>
      as.matrix()
    mat_part2 <- data$df_observed |>
      filter(X1 <= 0) |>
      select(time, contains("X")) |>
      as.matrix()

    lpme_results <- list()
    # fit LPME algorithm and calculate reconstructions
    lpme_results[[1]] <- lpme(
      mat_part1,
      d,
      verbose = FALSE,
      print_plots = FALSE
    )
    lpme_reconstructions_part1 <- calculate_lpme_reconstructions(
      lpme_results[[1]],
      mat_part1
    )

    lpme_results[[2]] <- lpme(
      mat_part2,
      d,
      verbose = FALSE,
      print_plots = FALSE
    )
    lpme_reconstructions_part2 <- calculate_lpme_reconstructions(
      lpme_results[[2]],
      mat_part2
    )

    lpme_reconstructions <- rbind(
      lpme_reconstructions_part1,
      lpme_reconstructions_part2
    )

    # fit PME and appropriate principal curve/surface algorithms and
    # calculate reconstructions for each
    pme_result_list <- list()
    pc_result_list <- list()

    pme_reconstruction_list <- list()
    pc_reconstruction_list <- list()

    smoothing_options <- c("smooth_spline", "lowess", "periodic_lowess")

    time_values <- unique(data$df$time)
    for (time_idx in seq_along(time_values)) {
      temp_data_part1 <- mat_part1[mat_part1[, 1] == time_values[time_idx], -1]
      temp_data_part2 <- mat_part2[mat_part2[, 1] == time_values[time_idx], -1]

      temp_pme_results <- list()
      temp_pme_results[[1]] <- pme(temp_data_part1, d = d)
      temp_pme_results[[2]] <- pme(temp_data_part2, d = d)
      pme_result_list[[time_idx]] <- temp_pme_results

      pme_reconstructions_part1 <- calculate_pme_reconstructions(
        pme_result_list[[time_idx]][[1]],
        temp_data_part1
      )
      pme_reconstructions_part2 <- calculate_pme_reconstructions(
        pme_result_list[[time_idx]][[2]],
        temp_data_part2
      )
      pme_reconstruction_list[[time_idx]] <- rbind(
        pme_reconstructions_part1,
        pme_reconstructions_part2
      )
      pme_reconstruction_list[[time_idx]] <- cbind(
        time_values[time_idx],
        pme_reconstruction_list[[time_idx]]
      )

      pc_result_list[[time_idx]] <- list()
      principal_surface_part1 <- prinSurf(temp_data_part1)
      surface_mse_part1 <- map(
        seq_along(principal_surface_part1),
        ~ principal_surface_part1[[.x]]$MSE
      ) |>
        reduce(c)
      opt_surface_part1 <- which.min(surface_mse_part1)
      pc_result_part1 <- principal_surface_part1[[opt_surface_part1 + 2]]
      pc_result_list[[time_idx]][[1]] <- pc_result_part1

      pc_reconstructions_part1 <- cbind(
        time_values[time_idx],
        pc_result_part1$PS
      )

      principal_surface_part2 <- prinSurf(temp_data_part2)
      surface_mse_part2 <- map(
        seq_along(principal_surface_part2),
        ~ principal_surface_part2[[.x]]$MSE
      ) |>
        reduce(c)
      opt_surface_part2 <- which.min(surface_mse_part2)
      pc_result_part2 <- principal_surface_part2[[opt_surface_part2 + 2]]
      pc_result_list[[time_idx]][[2]]

      pc_reconstructions_part2 <- cbind(
        time_values[time_idx],
        pc_result_part2$PS
      )

      pc_reconstruction_list[[time_idx]] <- rbind(
        pc_reconstructions_part1,
        pc_reconstructions_part2
      )
    }
  }

  pme_reconstructions <- reduce(pme_reconstruction_list, rbind)
  if (length(pc_reconstruction_list) > 0) {
    pc_reconstructions <- reduce(pc_reconstruction_list, rbind)
  } else {
    pc_reconstructions <- NULL
  }

  lpme_out <- list(
    lpme = lpme_results,
    reconstructions = lpme_reconstructions
  )

  pme_out <- list(
    pme = pme_result_list,
    reconstructions = pme_reconstructions
  )

  pc_out <- list(
    pc = pc_result_list,
    reconstructions = pc_reconstructions
  )

  model_out <- list(
    lpme = lpme_out,
    pme = pme_out,
    pc = pc_out
  )

  return(model_out)
}
