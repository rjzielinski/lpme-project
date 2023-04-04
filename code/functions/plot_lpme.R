plot_lpme <- function(x.obs, f, tnew, d_new, D_new, time_points) {
  time_vals <- seq(
    min(time_points),
    max(time_points),
    0.05
  )
  # pred_grid <- calc_tnew(centers, tnew, sol, I, d_new, gamma)
  # r_bounds <- colMinsMaxs(pred_grid)
  # r_list <- list()
  # for (idx in 1:dim(r_bounds)[2]) {
  #   r_list[[idx]] <- seq(
  #     r_bounds[1, idx],
  #     r_bounds[2, idx],
  #     length.out = sqrt(nrow(centers))
  #   )
  # }
  # r_mat <- as.matrix(expand.grid(r_list))
  # pred_grid <- r_mat
  new_pred <- matrix(ncol = ncol(tnew) + 1)
  for (time in time_vals) {
    new_pred <- rbind(new_pred, cbind(time, tnew))
  }
  pred_grid <- new_pred[-1, ]
  f_pred <- map(
    1:nrow(pred_grid),
    ~ f(pred_grid[.x, ])
  ) %>%
    reduce(rbind)

  f_pred_full <- cbind(pred_grid, f_pred)

  # if (D_new == 3) {
  if (D_new == 2) {
    plt <- plot_ly(
      x = f_pred_full[, d_new + 3],
      y = f_pred_full[, d_new + 4],
      # x = f_pred_full[, d_new + 1],
      # y = f_pred_full[, d_new + 2],
      z = f_pred_full[, 1],
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = 3
      )
    ) %>%
      add_markers(
        x = x.obs[, 2],
        y = x.obs[, 3],
        z = x.obs[, 1],
        opacity = 0.15
      )
    print(plt)
  } else {
    if (D_new >= 3) {
    # if (D_new >= 4) {
    plt <- plot_ly(
      x = f_pred_full[, d_new + 3],
      # x = f_pred_full[, d_new + 1],
      y = f_pred_full[, d_new + 4],
      # y = f_pred_full[, d_new + 2],
      z = f_pred_full[, d_new + 5],
      # z = f_pred_full[, d_new + 3],
      frame = f_pred_full[, d_new + 2],
      # frame = f_pred_full[, d_new],
      type = "scatter3d",
      mode = "markers",
      opacity = 1,
      marker = list(size = 3)
    ) %>%
      add_markers(
        x = x.obs[, 2],
        y = x.obs[, 3],
        z = x.obs[, 4],
        frame = x.obs[, 1],
        opacity = 0.2
      ) %>%
      layout(
        scene = list(
          xaxis = list(
            range = list(
              min(x.obs[, 2]),
              max(x.obs[, 2])
            )
          ),
          yaxis = list(
            range = list(
              min(x.obs[, 3]),
              max(x.obs[, 3])
            )
          ),
          zaxis = list(
            range = list(
              min(x.obs[, 4]),
              max(x.obs[, 4])
            )
          )
        )
      )
    print(plt)
    }
  }
}
