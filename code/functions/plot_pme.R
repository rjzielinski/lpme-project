plot_pme <- function(x.obs, centers, sol, tnew, I, d, lambda) {
  pred_grid <- calc_tnew(centers, tnew, sol, I, d, lambda)
  r_bounds <- colMinsMaxs(pred_grid)
  r_list <- list()
  for (idx in 1:dim(r_bounds)[2]) {
    r_list[[idx]] <- seq(
      r_bounds[1, idx],
      r_bounds[2, idx],
      length.out = nrow(centers)
    )
  }
  r_mat <- as.matrix(expand.grid(r_list))

  pred_grid <- r_mat
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

  if (dim(x.obs)[2] == 2) {
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
  } else if (dim(x.obs)[2] >= 3) {
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
