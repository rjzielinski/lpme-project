interior_identification <- function(embedding_map, coefs, params, x, reference) {
  orientation_c <- get_orientation(embedding_map, coefs, params, reference)[1]
  orientation_x <- map(
    1:nrow(x), 
    ~ get_orientation(embedding_map, coefs, params, unlist(x[.x, ]))
  ) %>% 
    unlist()
  return(orientation_x == orientation_c)
}

  get_orientation <- function(embedding_map, coefs, params, x) {
  d <- ncol(params)
  n <- nrow(params)

  augment_x <- c(x, cart2sph(x))
  x_params <- projection_pme(augment_x, embedding_map, params[1, ])
  x_proj <- embedding_map(x_params)

  
  param_diff <- -1 * sweep(params, 2, x_params)
  param_diff_norms <- wordspace::rowNorms(param_diff)

  jacobian <- t(coefs[1:n, ]) %*% (param_diff * ((2 * log(param_diff_norms)) + 1)) +
    t(coefs[(n + 1):(n + d + 1), ][-1, ])

  normal_vector <- c(
    (jacobian[2, 1] * jacobian[3, 2]) - (jacobian[3, 1] * jacobian[2, 2]),
    (jacobian[3, 1] * jacobian[1, 2]) - (jacobian[1, 1] * jacobian[3, 2]),
    (jacobian[1, 1] * jacobian[2, 2]) - (jacobian[2, 1] * jacobian[1, 2])
  ) %>%
    matrix(ncol = 1)

  orientation <- sign(t(matrix(x_proj - augment_x, nrow = 1)[, 1:3]) %*% normal_vector)
  return(orientation)
}
