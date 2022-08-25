library(tidyverse)
library(plotly)

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

time_vals <- seq(0, 2, 0.25)

df_list <- lapply(
  time_vals, 
  sim_D2d1_case1, 
  vertical_multiplier = 1, 
  horizontal_multiplier = 1, 
  noise = 0.15, 
  time_noise = 1
)

data_points <- reduce(df_list, rbind)

sim_result <- long_pme(data_points, 1)
