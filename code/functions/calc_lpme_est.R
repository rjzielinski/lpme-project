calc_lpme_est <- function(lpme, df) {
  clusters <- lpme$clusters
  opt_run <- which.min(lpme$MSD)
  tnew_init <- lpme$TNEW[[opt_run]][clusters, ]
  tnew <- map(
    1:nrow(df),
    ~ projection_lpme(df[.x, ], lpme$embedding_map, tnew_init[.x, ])
  ) %>%
    reduce(rbind)
  results <- map(
    1:nrow(df),
    ~ lpme$embedding_map(tnew[.x, ])
  ) %>%
    reduce(rbind)
  return(results)
}
