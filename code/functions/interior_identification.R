interior_identification <- function(embedding_map, coefs, params, x, reference) {
  orientation_x <- get_orientation(embedding_map, coefs, params, x)
  orientation_c <- get_orientation(embedding_map, coefs, params, reference)
  return(orientation_x == orientation_c)
}

get_orientation <- function(embedding_map, coefs, params, x) {
  # This function identifies the interior of a manifold with d = 2 and D = 3,
  # estimated using PME or LPME.
  # Args:
  #   embedding_map: the embedding function estimated by PME
  #   coefs: the coefficients of the embedding map
  #   params: the parameters of the embedding map
  #   x: the point to be identified
  # Returns: a logical value indicating whether x is in the interior of the manifold or not

  d <- ncol(params)
  n <- nrow(params)

  augment_x <- c(x, cart2sph(x))
  x_params <- projection_pme(augment_x, embedding_map, params[1, ])
  x_proj <- embedding_map(x_params)

  param_diff <- map(1:nrow(params), ~x_params - params[.x, ]) %>%
    reduce(rbind)
  # param_diff <- x_params - params
  # norm_param_diff <- norm_euclidean(param_diff)
  param_diff_norms <- map(
    1:nrow(param_diff),
    ~norm_euclidean(param_diff[.x, ])
  ) %>%
    reduce(c)

  jacobian <- t(coefs[1:n, ]) %*% (param_diff * ((2 * log(param_diff_norms)) + 1)) +
    t(coefs[(n + 1):(n + d + 1), ][-1, ])

  normal_vector <- c(
    (jacobian[2, 1] * jacobian[3, 2]) - (jacobian[3, 1] * jacobian[2, 2]),
    (jacobian[3, 1] * jacobian[1, 2]) - (jacobian[1, 1] * jacobian[3, 2]),
    (jacobian[1, 1] * jacobian[2, 2]) - (jacobian[2, 1] * jacobian[1, 2])
  ) %>%
    matrix(ncol = 1)

  # normal_vector <- pracma::crossn(t(jacobian))

  orientation <- sign(t(matrix(x_proj - augment_x, nrow = 1)[, 1:3]) %*% normal_vector)
  return(orientation)
}
