long_pme <- function(df, d, tuning.para.seq = exp(-15:5), alpha = 0.05, max.comp = 100, epsilon = 0.05, max.iter = 100) {
  # df is an N x (D + 1) matrix, with the first column corresponding
  # to the time point at which each observation was collected
  # this matrix should include the observations from all time points
  
  source("PME_recode.R")
  # source("Principal_Manifold_Estimation.R")
  
  time_points <- df[, 1] %>% 
    unique()
  
  pme_results <- list()
  funcs <- list()
  x_test <- list()
  
  for (idx in 1:length(time_points)) {
    df_temp <- df[df[, 1] == time_points[idx], ]
    pme_results[[idx]] <- PME(
      x.obs = df_temp[, -1],
      d = d
    )
    funcs[[idx]] <- pme_results[[idx]]$embedding.map
    t_test <- seq(from = -100, to = 100, by = 0.05)
    t_length <- length(t_test)
    x_test_temp <- map(t_test, ~ funcs[[idx]](.x)) %>% 
      unlist() %>% 
      matrix(nrow = t_length, byrow = TRUE)
    x_test[[idx]] <- cbind(time_points[idx], x_test_temp)
  }
  
  x_test <- reduce(x_test, rbind) 
  x_test_df <- data.frame(x_test)
  names(x_test_df) <- c("time", "x", "y")
  x_test_df <- x_test_df %>% 
    mutate(time = as.factor(time))
  
  # plot_ly(
  #   x_test_df,
  #   x = ~x,
  #   y = ~y,
  #   z = ~time,
  #   type = "scatter3d",
  #   mode = "lines"
  # )
  
  D_new <- dim(x_test)[2]
  d_new <- dim(x_test)[2] - 1
  n_new <- dim(x_test)[1]
  gamma <- 4 - d_new
  
  
}