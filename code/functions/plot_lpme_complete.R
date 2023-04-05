plot_lpme_complete <- function(df, result) {
  opt_result <- which.min(result$MSD)
  # t_new <- result$TNEW[[opt_result]]
  coefs <- result$coefs[[opt_result]]
  d_new <- ncol(coefs)
  D_new <- ncol(df) - 1
  time_points <- unique(df[,  1])
  plot_lpme(df, result$embedding_map, coefs, d_new, D_new, time_points)
}
