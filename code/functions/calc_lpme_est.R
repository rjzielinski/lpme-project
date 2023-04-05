calc_lpme_est <- function(result, df) {
  opt_idx <- which.min(result$MSD)
  x_vals <- result$x_vals[[opt_idx]]
  tnew_init <- result$TNEW[[opt_idx]]

  nearest_x <- map(
    1:nrow(df),
    ~ which.min(apply(x_vals[x_vals[, 1] == df[.x, 1], ], 1, dist_euclideanC, y = df[.x, ]))
  ) %>%
    reduce(c)

  init_param <- map(
    1:nrow(df),
    ~ tnew_init[tnew_init[, 1] == df[.x, 1], ][nearest_x[.x], ]
  ) %>%
    reduce(rbind)

  tnew <- map(
    1:nrow(df),
    ~ projection_lpme(df[.x, ], result$embedding_map, init_param[.x, ])
  ) %>%
    reduce(rbind)
  results <- map(
    1:nrow(df),
    ~ result$embedding_map(tnew[.x, ])
  ) %>%
    reduce(rbind)
  return(results)
}
