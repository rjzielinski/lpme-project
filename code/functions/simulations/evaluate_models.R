evaluate_models <- function(data, models, case, d, D) {
  require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  require(purrr, quietly = TRUE, warn.conflicts = FALSE)

  true_values <- data$df_true |>
    select(time, contains("X")) |>
    as.matrix()

  observed_data <- data$df_observed |>
    select(time, contains("X")) |>
    as.matrix()

  lpme_error <- map(
    seq_len(nrow(true_values)),
    ~ dist_euclidean(
      true_values[.x, ],
      models$lpme$reconstructions[.x, 1:(1 + D)]
    )^2
  ) |>
    unlist() |>
    mean()

  pme_error <- map(
    seq_len(nrow(true_values)),
    ~ dist_euclidean(
      true_values[.x, ], models$pme$reconstructions[.x, 1:(1 + D)]
    )^2
  ) |>
    unlist() |>
    mean()

  if (!is.null(models$pc$reconstructions)) {
    principal_curve_error <- map(
      seq_len(nrow(true_values)),
      ~ dist_euclidean(
        true_values[.x, ], models$pc$reconstructions[.x, 1:(1 + D)]
      )^2
    ) |>
      unlist() |>
      mean()
  } else {
    principal_curve_error <- NA
  }

  if (!is.null(models$lpme_part$reconstructions)) {
    lpme_part_error <- map(
      seq_len(nrow(true_values)),
      ~ dist_euclidean(
        true_values[.x, ], models$lpme_part$reconstructions[.x, 1:(1 + D)]
      )^2
    ) |>
      unlist() |>
      mean()
  } else {
    lpme_part_error <- NA
  }

  if (!is.null(models$pme_part$reconstructions)) {
    pme_part_error <- map(
      seq_len(nrow(true_values)),
      ~ dist_euclidean(
        true_values[.x, ], models$pme_part$reconstructions[.x, 1:(1 + D)]
      )^2
    ) |>
      unlist() |>
      mean()
  } else {
    pme_part_error <- NA
  }

  if (!is.null(models$pc_part$reconstructions)) {
    principal_curve_part_error <- map(
      seq_len(nrow(true_values)),
      ~ dist_euclidean(
        true_values[.x, ], models$pc_part$reconstructions[.x, 1:(1 + D)]
      )^2
    ) |>
      unlist() |>
      mean()
  } else {
    principal_curve_part_error <- NA
  }

  data_error <- map(
    seq_len(nrow(true_values)),
    ~ dist_euclidean(true_values[.x, ], observed_data[.x, ])^2
  ) |>
    unlist() |>
    mean()

  error_list <- list(
    lpme_error = lpme_error,
    pme_error = pme_error,
    principal_curve_error = principal_curve_error,
    lpme_part_error = lpme_part_error,
    pme_part_error = pme_part_error,
    principal_curve_part_error = principal_curve_part_error,
    data_error = data_error
  )

  return(error_list)
}
