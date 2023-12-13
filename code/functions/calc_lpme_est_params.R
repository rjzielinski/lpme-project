calc_lpme_est_params <- function(result, df) {
  opt_idx <- which.min(result$msd)
  r_full <- expand_grid(result$times, result$initial_parameterization) %>%
    as.matrix()
  # x_vals <- map(
  #   1:nrow(r_full),
  #   ~ embed(result, r_full[.x, ])
  # ) %>%
  #   reduce(rbind)
  x_vals <- map(
    1:nrow(r_full),
    ~ result$embedding_map(r_full[.x, ])
  ) %>%
    reduce(rbind)
  tnew_init <- result$optimal_parameterization

  nearest_x <- map(
    1:nrow(df),
    ~ which.min(apply(x_vals[x_vals[, 1] == df[.x, 1], ], 1, dist_euclidean, y = df[.x, ]))
  ) %>%
    reduce(c)

  init_param <- map(
    1:nrow(df),
    ~ tnew_init[tnew_init[, 1] == df[.x, 1], ][nearest_x[.x], ]
  ) %>%
    reduce(rbind)

  tnew <- map(
    1:nrow(df),
    ~ projection_lpme(df[.x, ], function(x) result$embedding_map(x), init_param[.x, ])
  ) %>%
    reduce(rbind)
  results <- map(
    1:nrow(df),
    ~ result$embedding_map(tnew[.x, ])
  ) %>%
    reduce(rbind)
  return(list(results = results, params = tnew))
}
