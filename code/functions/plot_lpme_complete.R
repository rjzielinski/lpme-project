plot_lpme_complete <- function(df, result) {
  opt_result <- which.min(result$MSD)
  t_new <- result$TNEW[[opt_result]]
  d_new <- ncol(t_new)
  D_new <- ncol(df)
  time_points <- unique(df[,  1])
  plot_lpme(df, result$embedding_map, t_new, d_new, D_new, time_points)
}