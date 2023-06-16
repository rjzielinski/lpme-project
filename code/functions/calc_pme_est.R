calc_pme_est <- function(pme, df) {
  clusters <- pme$knots$cluster
  opt_run <- which.min(pme$MSD)
  tnew_init <- pme$parameterization[[opt_run]][clusters, ] %>%
    matrix(nrow = nrow(df), byrow = TRUE)
  sol_init <- pme$coefs[[opt_run]][clusters, ] %>%
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
