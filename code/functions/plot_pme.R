plot_pme <- function(time_points, x.obs, sol, tnew, I, d, lambda) {
  time_vals <- seq(
    min(time_points),
    max(time_points),
    0.1
  )

  r_vals <- seq(
    from = -10,
    to = 10,
    by = 1
  )
  r_list <- lapply(numeric(d - 1), function(x) r_vals)
  r_mat <- as.matrix(expand.grid(r_list))
  if (nrow(r_mat) > 0) {
    pred_grid <- expand_grid(r_vals, r_mat) %>%
      as.matrix()
  } else {
    pred_grid <- matrix(r_vals)
  }

  f_pred <- map(
    1:nrow(pred_grid),
    ~ fNew(
      unlist(as.vector(pred_grid[.x, ])),
      sol,
      tnew,
      I,
      d,
      lambda
    ) %>%
      as.vector()
  ) %>%
    unlist() %>%
    matrix(nrow = nrow(pred_grid), byrow = TRUE)

  idx_inrange <- matrix(nrow = dim(f_pred)[1], ncol = dim(f_pred)[2])
  for (dim_idx in 1:dim(f_pred)[2]) {
    idx_range <- max(x.obs[, dim_idx]) - min(x.obs[, dim_idx])
    idx_min <- min(x.obs[, dim_idx]) - (0.2 * idx_range)
    idx_max <- max(x.obs[, dim_idx]) + (0.2 * idx_range)
    idx_inrange[, dim_idx] <- (f_pred[, dim_idx] > idx_min) &
      (f_pred[, dim_idx] < idx_max)
  }

  r_inrange <- rowSums(idx_inrange) == dim(f_pred)[2]
  if (sum(r_inrange) == 0) {
    r_inrange <- rowSums(idx_inrange) > 0
  }
  if (sum(r_inrange) == 0) {
    r_min <- -10
    r_max <- 10
  } else if (dim(pred_grid)[2] > 1) {
    r_min <- min(pred_grid[r_inrange, -1])
    r_min <- ifelse(is.na(r_min) | is.infinite(r_min), -10, r_min)
    r_max <- max(pred_grid[r_inrange, -1])
    r_max <- ifelse(is.na(r_max) | is.infinite(r_max), 10, r_max)
  } else {
    r_min <- min(pred_grid[r_inrange, ])
    r_min <- ifelse(is.na(r_min) | is.infinite(r_min), -10, r_min)
    r_max <- max(pred_grid[r_inrange, ])
    r_max <- ifelse(is.na(r_max) | is.infinite(r_max), 10, r_max)
  }

  r_vals <- seq(
    r_min,
    r_max,
    by = (r_max - r_min) / 40
  )
  r_list <- lapply(numeric(d - 1), function(x) r_vals)
  r_mat <- as.matrix(expand.grid(r_list))
  if (nrow(r_mat) > 0) {
    pred_grid <- expand_grid(r_vals, r_mat) %>%
      as.matrix()
  } else {
    pred_grid <- as.matrix(r_vals)
  }
  f_pred <- map(
    1:nrow(pred_grid),
    ~ fNew(
      unlist(as.vector(pred_grid[.x, ])),
      sol,
      tnew,
      I,
      d,
      lambda
    ) %>%
      as.vector()
  ) %>%
    unlist() %>%
    matrix(nrow = nrow(pred_grid), byrow = TRUE)

  f_pred_full <- cbind(pred_grid, f_pred)

  if (D == 2) {
    plt <- ggplot() +
      geom_point(
        aes(
          x =x.obs[, 1],
          y = x.obs[, 2]
        ),
        alpha = 0.5
      ) +
      geom_point(
        aes(
          x = f_pred_full[, d + 1],
          y = f_pred_full[, d + 2]
        ),
        color = "red"
      )
    print(plt)
  } else if (D == 3) {
    plt <- plot_ly(
      x = f_pred_full[, d + 1],
      y = f_pred_full[, d + 2],
      z = f_pred_full[, d + 3],
      type = "scatter3d",
      mode = "markers",
      opacity = 0.5
    ) %>%
      add_markers(
        x = x.obs[, 1],
        y = x.obs[, 2],
        z = x.obs[, 3],
        opacity = 0.15
      )
    print(plt)
  }
}