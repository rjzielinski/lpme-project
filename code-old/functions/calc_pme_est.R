calc_pme_est <- function(pme, df) {
  center_order <- order(pme$knots$centers[, 1])
  centers <- pme$knots$centers[center_order, ]

  nearest_clusters <- map(
    1:nrow(df),
    ~ which.min(as.vector(apply(centers, 1, dist_euclidean, y = df[.x, ])))
  ) %>%
    reduce(c)

  opt_run <- which.min(pme$MSD)
  tnew_init <- pme$parameterization[[opt_run]][nearest_clusters, ] %>%
    matrix(nrow = nrow(df), byrow = TRUE)

  tnew <- map(
    1:nrow(df),
    ~ projection_pme(df[.x, ], pme$embedding_map, tnew_init[.x, ])
  ) %>%
    reduce(rbind)
  results <- map(
    1:nrow(df),
    ~ pme$embedding_map(tnew[.x, ])
  ) %>%
    reduce(rbind)
  return(results)
}
