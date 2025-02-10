interior_identification <- function(embedding_map, coefs, params, x, ref) {
  require(pme, quietly = TRUE, warn.conflicts = FALSE)
  require(pracma, quietly = TRUE, warn.conflicts = FALSE)
  require(purrr, quietly = TRUE, warn.conflicts = FALSE)
  require(wordspace, quietly = TRUE, warn.conflicts = FALSE)

  param_centers <- map(
    seq_len(nrow(params)),
    ~ embedding_map(params[.x, ])
  ) |>
    reduce(rbind)

  orientation_c <- get_orientation(
    embedding_map,
    coefs,
    param_centers,
    params,
    ref
  )[1]
  orientation_x <- map(
    seq_len(nrow(x)),
    ~ get_orientation(
      embedding_map,
      coefs,
      param_centers,
      params,
      unlist(x[.x, ])
    )
  ) %>%
    unlist()
  return(orientation_x == orientation_c)
}

get_orientation <- function(embedding_map, coefs, centers, params, x) {
  d <- ncol(params)
  n <- nrow(params)

  augment_x <- c(x, cart2sph(x))[-6]
  nearest_center <- map(
    seq_len(nrow(centers)),
    ~ dist_euclidean(augment_x, as.vector(centers[.x, ]))
  ) |>
    reduce(c) |>
    which.min()
  x_params <- projection_pme(augment_x, embedding_map, params[nearest_center, ])
  x_proj <- embedding_map(x_params)


  param_diff <- -1 * sweep(params, 2, x_params)
  param_diff_norms <- rowNorms(param_diff)

  jacobian <- t(coefs[1:n, ]) %*%
    (param_diff * ((2 * log(param_diff_norms)) + 1)) +
    t(coefs[(n + 1):(n + d + 1), ][-1, ])

  normal_vector <- cross(
    as.vector(jacobian[1:3, 1]),
    as.vector(jacobian[1:3, 2])
  ) %>%
    matrix(ncol = 1)

  orientation <- sign(
    t(matrix(x_proj - augment_x, nrow = 1)[, 1:3]) %*% normal_vector
  )
  return(orientation)
}
