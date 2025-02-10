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

  augment_x <- cbind(x, cart2sph(x))[, -6]
  augment_c <- c(ref, cart2sph(ref))[-6]
  nearest_centers <- vector(mode = "numeric", length = nrow(augment_x))
  for (row_idx in seq_len(nrow(augment_x))) {
    nearest_centers[row_idx] <- map(
      seq_len(nrow(param_centers)),
      ~ dist_euclidean(
        as.vector(augment_x[row_idx, ]),
        as.vector(param_centers[.x, ])
      )
    ) |>
      reduce(c) |>
      which.min()
  }

  nearest_center_c <- map(
    seq_len(nrow(param_centers)),
    ~ dist_euclidean(
      as.vector(augment_c),
      as.vector(param_centers[.x, ])
    )
  ) |>
    reduce(c) |>
    which.min()

  c_param <- projection_pme(
    augment_c,
    embedding_map,
    params[nearest_center_c, ]
  )
  c_proj <- embedding_map(c_param)
  c_diff <- augment_c - c_proj

  x_params <- map(
    seq_len(nrow(augment_x)),
    ~ projection_pme(
      as.vector(augment_x[.x, ]),
      embedding_map,
      params[nearest_centers[.x], ]
    )
  ) |>
    reduce(rbind)
  x_proj <- map(
    seq_len(nrow(x_params)),
    ~ embedding_map(as.vector(x_params[.x, ]))
  ) |>
    reduce(rbind)
  x_diff <- augment_x - x_proj

  angles <- map(
    seq_len(nrow(x_diff)),
    ~ acos(
      (x_diff[.x, ] %*% c_diff)[1] /
        (norm_euclidean(x_diff[.x, ]) * norm_euclidean(c_diff))
    )
  ) |>
    reduce(c)

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
