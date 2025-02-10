preprocess_data <- function(df, case, d, D) {
  require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  require(pracma, quietly = TRUE, warn.conflicts = FALSE)


  ##### PREPROCESSING *****

  # standardize observations to fall between -1 and 1 in each dimension
  # use maximum observed values to standardize both observations and true values
  time_values <- unique(df$time)

  for (dim_idx in 1:D) {
    obs_dim_idx <- 1 + d + dim_idx
    true_dim_idx <- 1 + d + D + dim_idx
    col_max <- max(abs(df[obs_dim_idx]))
    df[, obs_dim_idx] <- df[, obs_dim_idx] / col_max
    df[, true_dim_idx] <- df[, true_dim_idx] / col_max
  }
  df$time <- df$time / max(df$time)
  time_values <- time_values / max(time_values)

  # where, if anywhere, should case 7 (swiss roll) be included?
  polar_cases <- c(3)
  spherical_cases <- c(8)

  if (case %in% polar_cases) {
    observed_polar <- df |>
      select(contains("X")) |>
      select(-contains("true")) |>
      as.matrix() |>
      cart2pol()

    true_polar <- df |>
      select(contains("true")) |>
      as.matrix() |>
      cart2pol()

    observed_polar <- observed_polar[, 1]
    true_polar <- true_polar[, 1]

    df$phi <- observed_polar
    df$phi_true <- true_polar
  } else if (case %in% spherical_cases) {
    observed_spherical <- df |>
      select(contains("X")) |>
      select(-contains("true")) |>
      as.matrix() |>
      cart2sph()

    true_spherical <- df |>
      select(contains("true")) |>
      as.matrix() |>
      cart2sph()

    observed_spherical_df <- observed_spherical[, -3] |>
      as.data.frame()
    names(observed_spherical_df) <- c("theta", "phi")

    true_spherical_df <- true_spherical[, -3] |>
      as.data.frame()
    names(true_spherical_df) <- c("theta_true", "phi_true")

    df <- df |>
      bind_cols(observed_spherical_df) |>
      bind_cols(true_spherical_df)
  }

  df_observed <- df |>
    select(time, contains("X"), contains("phi"), contains("theta")) |>
    select(-contains("true"))

  df_true <- df |>
    select(time, contains("true"))

  data <- list(
    df = df,
    df_observed = df_observed,
    df_true = df_true
  )
  return(data)
}
