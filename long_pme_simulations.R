library(tidyverse)
library(plotly)
library(profvis)

source("long_pme.R")

sim_D2d1_case1 <- function(time_val, vertical_multiplier, horizontal_multiplier, noise, time_noise) { 
  I <- 1000
  t <- rnorm(I, mean = 0, sd = 1)
  horizontal_noise <- rnorm(1, mean = 0, sd = time_noise)
  vertical_noise <- rnorm(1, mean = 0, sd = time_noise)
  sd.noise <- noise
  e1 <- rnorm(I, mean = 0, sd = sd.noise)
  e2 <- rnorm(I, mean = 0, sd = sd.noise)
  X <- matrix(NA, nrow = I, ncol = 2)
  manifold <- function(tau, time_val, vertical_multiplier, horizontal_multiplier, vertical_noise, horizontal_noise) {
    return(
      c(
        tau + (horizontal_multiplier * sin(time_val)),
        # tau + horizontal_noise,
        sin(tau + pi / 2) + (vertical_multiplier * sin(time_val))
        # sin(tau + pi / 2) + vertical_noise
      )
    )
  }
  
  X <- map(
    t, 
    ~ manifold(
      .x, 
      time_val, 
      vertical_multiplier, 
      horizontal_multiplier, 
      vertical_noise, 
      horizontal_noise
    )
  ) %>% 
    unlist() %>% 
    matrix(ncol = 2, byrow = TRUE)
  data.points <- X + cbind(e1, e2)
  data.points <- cbind(time_val, data.points)
  return(data.points)
}

sim_D2d1_case2 <- function(time_val, vertical_multiplier, horizontal_multiplier, noise, time_noise) { 
  I <- 1000
  t <- rnorm(I, mean = 0, sd = 1)
  horizontal_noise <- rnorm(1, mean = 0, sd = time_noise)
  vertical_noise <- rnorm(1, mean = 0, sd = time_noise)
  sd.noise <- noise
  e1 <- rnorm(I, mean = 0, sd = sd.noise)
  e2 <- rnorm(I, mean = 0, sd = sd.noise)
  X <- matrix(NA, nrow = I, ncol = 2)
  manifold <- function(tau, time_val, vertical_multiplier, horizontal_multiplier, vertical_noise, horizontal_noise) {
    return(
      c(
        # cos(tau) + (horizontal_multiplier * sin(time_val)),
        cos(tau) + horizontal_noise,
        # sin(tau) + (vertical_multiplier * sin(time_val))
        sin(tau) + vertical_noise
      )
    )
  }
  
  X <- map(
    t, 
    ~ manifold(
      .x, 
      time_val, 
      vertical_multiplier, 
      horizontal_multiplier, 
      vertical_noise, 
      horizontal_noise
    )
  ) %>% 
    unlist() %>% 
    matrix(ncol = 2, byrow = TRUE)
  data.points <- X + cbind(e1, e2)
  data.points <- cbind(time_val, data.points)
  return(data.points)
}

time_vals <- 0:3

set.seed(100)
df_list <- lapply(
  time_vals, 
  sim_D2d1_case1, 
  vertical_multiplier = 1, 
  horizontal_multiplier = 1, 
  noise = 0.15, 
  time_noise = 1
)

data_points <- reduce(df_list, rbind)

profvis({
  sim_result <- long_pme(data_points, 1)
})

time_vals <- seq(0, 3, 0.1)
r_vals <- seq(-10, 10, 0.1)
grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) & 
    (sim_pred[, dim_idx] < idx_max)
}
    
r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2][r_inrange, 1]))
r_max <- max(unlist(grid_mat[, 2][r_inrange, 1]))
r_vals <- seq(
  r_min,
  r_max,
  0.1
)

grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r", "x", "y")

plot_ly(
  sim_pred_full_df,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "markers"
)
