library(tidyverse)
library(plotly)

source("lpme.R")

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

# time_vals <- seq(0, 10, 0.25)
time_vals <- 0:5

df_list <- lapply(
  time_vals, 
  sim_D2d1_case1, 
  vertical_multiplier = 1, 
  horizontal_multiplier = 1, 
  noise = 0.15, 
  time_noise = 0.1
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 1)

df_list <- lapply(
  time_vals,
  sim_D2d1_case2,
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  noise = 0.15,
  time_noise = 0.1
)

data_points2 <- reduce(df_list, rbind)
sim_result2 <- long_pme(data_points2, 1)

time_test <- seq(0, 10, 0.1)
r_test <- seq(
  from = -10,
  to = 10,
  0.1
)

pars_test <- expand_grid(time_test, r_test)

sim_result1_pred <- apply(pars_test, 1, sim_result$embedding_map) %>% 
  matrix(ncol = 2, byrow = TRUE)
sim_result1_pred_df <- cbind(pars_test[, 1], sim_result1_pred) %>% 
  as_tibble()
names(sim_result1_pred_df) <- c("time", "x", "y")
sim_result1_pred_df <- sim_result1_pred_df %>% 
  mutate(
    idx = row_number(),
    in_range = (x > -7.5) & (x < 7.5) & (y > -5) & (y < 5)
  ) %>% 
  filter(in_range)

plot_ly(
  sim_result1_pred_df,
  x = ~y,
  y = ~x,
  z = ~time,
  type = "scatter3d"
)
