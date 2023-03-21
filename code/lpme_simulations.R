library(tidyverse)
library(plotly)
library(pracma)
library(profvis)
library(foreach)
library(doParallel)
library(doSNOW)
library(furrr)
library(progress)

source("code/functions/sim_data.R")
source("code/functions/calc_pme_est.R")
source("code/functions/calc_lpme_est.R")
source("code/pme.R")
source("code/lpme.R")

### Simulation Case 1

sim_error_case1 <- function(max_time, interval, noise, run = 1, print_plots = FALSE) {
  require(tidyverse)
  source("code/functions/sim_data.R")
  source("code/pme.R")
  source("code/lpme.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  time_vals <- seq(0, max_time, interval)
  sim_df <- lapply(
    time_vals,
    sim_data,
    case = 1,
    noise = 0.15,
    shape_noise = noise
  ) %>%
    reduce(rbind)

  lpme_result <- lpme(sim_df, 1, print_plots = print_plots, verbose = "MSD")
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  lpme_vals[, 1] <- sim_df[, 1]
  pme_result <- list()
  pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
  }
  pme_vals <- reduce(pme_vals, rbind)

  tau <- sim_df[, 2]
  true_vals <- matrix(nrow = nrow(sim_df), ncol = ncol(sim_df))
  true_vals[, 1] <- sim_df[, 1]
  true_vals[, 2] <- tau
  true_vals[, 3] <- sin(tau + (pi / 2))

  p <- plot_ly(
    x = lpme_vals[, 2],
    y = lpme_vals[, 3],
    z = lpme_vals[, 1],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 1.5)
  ) %>%
    add_markers(
      x = pme_vals[, 2],
      y = pme_vals[, 3],
      z = pme_vals[, 1]
    ) %>%
    add_markers(
      x = true_vals[, 2],
      y = true_vals[, 3],
      z = true_vals[, 1]
    )

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], pme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], lpme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  sim_case1 <- list(
    df = sim_df,
    times = time_vals,
    noise = noise,
    lpme_result = lpme_result,
    pme_results = pme_result,
    lpme_error = lpme_error,
    pme_error = pme_error
  )
  sim_dir <- "simulations/case1/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_noise_",
    str_pad(as.character(100 * noise), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case1,
    paste0(sim_dir, filename)
  )
  return(sim_case1)
}

### Simulation Case 2

sim_error_case2 <- function(max_time, interval, noise, run = 1, print_plots = FALSE) {
  require(tidyverse)
  source("code/functions/sim_data.R")
  source("code/pme.R")
  source("code/lpme.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  time_vals <- seq(0, max_time, interval)
  sim_df <- lapply(
    time_vals,
    sim_data,
    case = 2,
    noise = 0.15,
    shape_noise = noise
  ) %>%
    reduce(rbind)

  lpme_result <- lpme(sim_df, 1, print_plots = print_plots, verbose = "MSD")
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  lpme_vals[, 1] <- sim_df[, 1]
  pme_result <- list()
  pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
  }
  pme_vals <- reduce(pme_vals, rbind)

  tau <- sim_df[, 2]
  true_vals <- matrix(nrow = nrow(sim_df), ncol = ncol(sim_df))
  true_vals[, 1] <- sim_df[, 1]
  true_vals[, 2] <- tau
  true_vals[, 3] <- sin(tau)

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], pme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], lpme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  sim_case2 <- list(
    df = sim_df,
    times = time_vals,
    noise = noise,
    lpme_result = lpme_result,
    pme_results = pme_result,
    lpme_error = lpme_error,
    pme_error = pme_error
  )
  sim_dir <- "simulations/case2/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_noise_",
    str_pad(as.character(100 * noise), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case2,
    paste0(sim_dir, filename)
  )
  return(sim_case2)
}

### Simulation Case 3

sim_error_case3 <- function(max_time, interval, noise, run = 1, print_plots = FALSE) {
  require(tidyverse)
  source("code/functions/sim_data.R")
  source("code/pme.R")
  source("code/lpme.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  time_vals <- seq(0, max_time, interval)
  sim_df <- lapply(
    time_vals,
    sim_data,
    case = 3,
    noise = 0.15,
    shape_noise = noise
  ) %>%
    reduce(rbind)

  lpme_result <- lpme(sim_df, 1, print_plots = print_plots, verbose = "MSD")
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  lpme_vals[, 1] <- sim_df[, 1]
  pme_result <- list()
  pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
  }
  pme_vals <- reduce(pme_vals, rbind)

  tau <- acos(sim_df[, 2])
  true_vals <- matrix(nrow = nrow(sim_df), ncol = ncol(sim_df))
  true_vals[, 1] <- sim_df[, 1]
  true_vals[, 2] <- cos(tau)
  true_vals[, 3] <- sin(tau)

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], pme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], lpme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  sim_case3 <- list(
    df = sim_df,
    times = time_vals,
    noise = noise,
    lpme_result = lpme_result,
    pme_results = pme_result,
    lpme_error = lpme_error,
    pme_error = pme_error
  )
  sim_dir <- "simulations/case3/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_noise_",
    str_pad(as.character(100 * noise), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case3,
    paste0(sim_dir, filename)
  )
  return(sim_case3)
}

### Simulation Case 4

sim_error_case4 <- function(max_time, interval, noise, run = 1, print_plots = FALSE) {
  require(tidyverse)
  source("code/functions/sim_data.R")
  source("code/pme.R")
  source("code/lpme.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  time_vals <- seq(0, max_time, interval)
  sim_df <- lapply(
    time_vals,
    sim_data,
    case = 4,
    noise = 0.15,
    shape_noise = noise
  ) %>%
    reduce(rbind)

  lpme_result <- lpme(sim_df, 1, print_plots = print_plots, verbose = "MSD")
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  lpme_vals[, 1] <- sim_df[, 1]
  pme_result <- list()
  pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
  }
  pme_vals <- reduce(pme_vals, rbind)

  tau <- acos(sim_df[, 2])
  true_vals <- matrix(nrow = nrow(sim_df), ncol = ncol(sim_df))
  true_vals[, 1] <- sim_df[, 1]
  true_vals[, 2] <- cos(tau)
  true_vals[, 3] <- sin(tau)

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], pme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], lpme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  sim_case4 <- list(
    df = sim_df,
    times = time_vals,
    noise = noise,
    lpme_result = lpme_result,
    pme_results = pme_result,
    lpme_error = lpme_error,
    pme_error = pme_error
  )
  sim_dir <- "simulations/case4/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_noise_",
    str_pad(as.character(100 * noise), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case4,
    paste0(sim_dir, filename)
  )
  return(sim_case4)
}

### Simulation Case 5

sim_error_case5 <- function(max_time, interval, noise, run = 1, print_plots = FALSE) {
  require(tidyverse)
  source("code/functions/sim_data.R")
  source("code/pme.R")
  source("code/lpme.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  time_vals <- seq(0, max_time, interval)
  sim_df <- lapply(
    time_vals,
    sim_data,
    case = 5,
    noise = 0.15,
    shape_noise = noise
  ) %>%
    reduce(rbind)

  lpme_result <- lpme(sim_df, 1, print_plots = print_plots, verbose = "MSD")
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  lpme_vals[, 1] <- sim_df[, 1]
  pme_result <- list()
  pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
  }
  pme_vals <- reduce(pme_vals, rbind)

  tau <- sim_df[, 2]
  true_vals <- matrix(nrow = nrow(sim_df), ncol = ncol(sim_df))
  true_vals[, 1] <- sim_df[, 1]
  true_vals[, 2] <- tau
  true_vals[, 3] <- tau ^ 2
  true_vals[, 4] <- tau ^ 3

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], pme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], lpme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  sim_case5 <- list(
    df = sim_df,
    times = time_vals,
    noise = noise,
    lpme_result = lpme_result,
    pme_results = pme_result,
    lpme_error = lpme_error,
    pme_error = pme_error
  )
  sim_dir <- "simulations/case5/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_noise_",
    str_pad(as.character(100 * noise), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case5,
    paste0(sim_dir, filename)
  )
  return(sim_case5)
}

### Simulation Case 6

sim_error_case6 <- function(max_time, interval, noise, run = 1, print_plots = FALSE) {
  require(tidyverse)
  source("code/functions/sim_data.R")
  source("code/pme.R")
  source("code/lpme.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  time_vals <- seq(0, max_time, interval)
  sim_df <- lapply(
    time_vals,
    sim_data,
    case = 6,
    noise = 0.15,
    shape_noise = noise
  ) %>%
    reduce(rbind)

  lpme_result <- lpme(sim_df, 1, print_plots = print_plots, verbose = "MSD")
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  lpme_vals[, 1] <- sim_df[, 1]
  pme_result <- list()
  pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
  }
  pme_vals <- reduce(pme_vals, rbind)

  tau <- sim_df[, 2]
  true_vals <- matrix(nrow = nrow(sim_df), ncol = ncol(sim_df))
  true_vals[, 1] <- sim_df[, 1]
  true_vals[, 2] <- tau
  true_vals[, 3] <- cos(tau)
  true_vals[, 4] <- sin(tau)

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], pme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], lpme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  sim_case6 <- list(
    df = sim_df,
    times = time_vals,
    noise = noise,
    lpme_result = lpme_result,
    pme_results = pme_result,
    lpme_error = lpme_error,
    pme_error = pme_error
  )
  sim_dir <- "simulations/case6/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 2, side = "left", pad = "0"),
    "_noise_",
    str_pad(as.character(100 * noise), 2, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case6,
    paste0(sim_dir, filename)
  )
  return(sim_case6)
}

### Simulation Case 7

sim_error_case7 <- function(max_time, interval, noise, run = 1, print_plots = FALSE) {
  require(tidyverse)
  source("code/functions/sim_data.R")
  source("code/pme.R")
  source("code/lpme.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  time_vals <- seq(0, max_time, interval)
  sim_df <- lapply(
    time_vals,
    sim_data,
    case = 7,
    noise = 0.15,
    shape_noise = noise
  ) %>%
    reduce(rbind)

  lpme_result <- lpme(sim_df, 1, print_plots = print_plots, verbose = "MSD")
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  lpme_vals[, 1] <- sim_df[, 1]
  pme_result <- list()
  pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
  }
  pme_vals <- reduce(pme_vals, rbind)

  tau <- matrix(ncol = 2)
  tau[, 1] <- sim_df[, 2]
  tau[, 2] <- sim_df[, 3]
  true_vals <- matrix(nrow = nrow(sim_df), ncol = ncol(sim_df))
  true_vals[, 1] <- sim_df[, 1]
  true_vals[, 2] <- tau[, 1]
  true_vals[, 3] <- tau[, 2]
  true_vals[, 4] <- norm_euclidean(tau)^2

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], pme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], lpme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  sim_case7 <- list(
    df = sim_df,
    times = time_vals,
    noise = noise,
    lpme_result = lpme_result,
    pme_results = pme_result,
    lpme_error = lpme_error,
    pme_error = pme_error
  )
  sim_dir <- "simulations/case7/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_noise_",
    str_pad(as.character(100 * noise), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case7,
    paste0(sim_dir, filename)
  )
  return(sim_case7)
}

### Simulation Case 8

sim_error_case8 <- function(max_time, interval, noise, run = 1, print_plots = FALSE) {
  require(tidyverse)
  source("code/functions/sim_data.R")
  source("code/pme.R")
  source("code/lpme.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  time_vals <- seq(0, max_time, interval)
  sim_df <- lapply(
    time_vals,
    sim_data,
    case = 8,
    noise = 0.15,
    shape_noise = noise
  ) %>%
    reduce(rbind)

  lpme_result <- lpme(sim_df, 1, print_plots = print_plots, verbose = "MSD")
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  lpme_vals[, 1] <- sim_df[, 1]
  pme_result <- list()
  pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
  }
  pme_vals <- reduce(pme_vals, rbind)

  tau <- matrix(ncol = 2)
  tau[, 1] <- sqrt(sim_df[, 2] ^ 2 + sim_df[, 3] ^ 2)
  tau[, 2] <- sim_df[, 4]
  true_vals <- matrix(nrow = nrow(sim_df), ncol = ncol(sim_df))
  true_vals[, 1] <- sim_df[, 1]
  true_vals[, 2] <- tau[, 1] * cos(tau[, 1])
  true_vals[, 3] <- tau[, 1] * sin(tau[, 1])
  true_vals[, 4] <- tau[, 2]

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], pme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], lpme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  sim_case8 <- list(
    df = sim_df,
    times = time_vals,
    noise = noise,
    lpme_result = lpme_result,
    pme_results = pme_result,
    lpme_error = lpme_error,
    pme_error = pme_error
  )
  sim_dir <- "simulations/case8/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_noise_",
    str_pad(as.character(100 * noise), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case8,
    paste0(sim_dir, filename)
  )
  return(sim_case8)
}

### Simulation Case 9

sim_error_case9 <- function(max_time, interval, noise, run = 1, print_plots = FALSE) {
  require(tidyverse)
  source("code/functions/sim_data.R")
  source("code/pme.R")
  source("code/lpme.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  time_vals <- seq(0, max_time, interval)
  sim_df <- lapply(
    time_vals,
    sim_data,
    case = 9,
    noise = 0.15,
    shape_noise = noise
  ) %>%
    reduce(rbind)

  lpme_result <- lpme(sim_df, 1, print_plots = print_plots, verbose = "MSD")
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  lpme_vals[, 1] <- sim_df[, 1]
  pme_result <- list()
  pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
  }
  pme_vals <- reduce(pme_vals, rbind)

  tau <- matrix(ncol = 3)
  tau[, 1] <- sqrt(sim_df[, 2] ^ 2 + sim_df[, 3] ^ 2 + sim_df[, 4])
  tau[, 2] <- atan(sim_df[, 3] / sim_df[, 2])
  tau[, 3] <- acos(sim_df[, 4] / tau[, 1])
  true_vals <- matrix(nrow = nrow(sim_df), ncol = ncol(sim_df))
  true_vals[, 1] <- sim_df[, 1]
  true_vals[, 2] <- tau[, 1] * cos(tau[, 2]) * sin(tau[, 3])
  true_vals[, 3] <- tau[, 1] * sin(tau[, 2]) * sin(tau[, 3])
  true_vals[, 4] <- tau[, 1] * cos(tau[, 3])

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], pme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], lpme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  sim_case9 <- list(
    df = sim_df,
    times = time_vals,
    noise = noise,
    lpme_result = lpme_result,
    pme_results = pme_result,
    lpme_error = lpme_error,
    pme_error = pme_error
  )
  sim_dir <- "simulations/case9/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_noise_",
    str_pad(as.character(100 * noise), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case9,
    paste0(sim_dir, filename)
  )
  return(sim_case9)
}

### Simulation Case 10

sim_error_case10 <- function(max_time, interval, noise, run = 1, print_plots = FALSE) {
  require(tidyverse)
  source("code/functions/sim_data.R")
  source("code/pme.R")
  source("code/lpme.R")
  source("code/functions/calc_pme_est.R")
  source("code/functions/calc_lpme_est.R")
  time_vals <- seq(0, max_time, interval)
  sim_df <- lapply(
    time_vals,
    sim_data,
    case = 10,
    noise = 0.15,
    shape_noise = noise
  ) %>%
    reduce(rbind)

  lpme_result <- lpme(sim_df, 1, print_plots = print_plots, verbose = "MSD")
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  lpme_vals[, 1] <- sim_df[, 1]
  pme_result <- list()
  pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1, verbose = "none")
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
  }
  pme_vals <- reduce(pme_vals, rbind)

  tau <- matrix(ncol = 3)
  tau[, 1] <- sqrt(sim_df[, 2] ^ 2 + sim_df[, 3] ^ 2 + sim_df[, 4])
  tau[, 2] <- atan(sim_df[, 3] / sim_df[, 2])
  tau[, 3] <- acos(sim_df[, 4] / tau[, 1])
  true_vals <- matrix(nrow = nrow(sim_df), ncol = ncol(sim_df))
  true_vals[, 1] <- sim_df[, 1]
  true_vals[, 2] <- tau[, 1] * cos(tau[, 2]) * sin(tau[, 3])
  true_vals[, 3] <- tau[, 1] * sin(tau[, 2]) * sin(tau[, 3])
  true_vals[, 4] <- tau[, 1] * cos(tau[, 3])

  pme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], pme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  lpme_error <- map(
    1:nrow(true_vals),
    ~ dist_euclideanC(true_vals[.x, ], lpme_vals[.x, ])
  ) %>%
    unlist() %>%
    mean()

  sim_case10 <- list(
    df = sim_df,
    times = time_vals,
    noise = noise,
    lpme_result = lpme_result,
    pme_results = pme_result,
    lpme_error = lpme_error,
    pme_error = pme_error
  )
  sim_dir <- "simulations/case10/"
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir)
  }
  filename <- paste0(
    "duration_",
    str_pad(as.character(max_time), 2, side = "left", pad = "0"),
    "_interval_",
    str_pad(as.character(100 * interval), 3, side = "left", pad = "0"),
    "_noise_",
    str_pad(as.character(100 * noise), 3, side = "left", pad = "0"),
    "_run_",
    str_pad(as.character(run), 2, side = "left", pad = "0"),
    ".rds"
  )
  saveRDS(
    sim_case10,
    paste0(sim_dir, filename)
  )
  return(sim_case10)
}

noise_vals <- seq(0, 2, 0.25)
max_times <- c(1, 2, 5, 8, 10)
intervals <- c(0.1, 0.25, 0.5, 1)
replicates <- 1:4

param_grid <- expand.grid(noise_vals, max_times, intervals, replicates)

plan(multisession, workers = availableCores() / 2)
# plan(sequential)
set.seed(1286)
pb <- progress_bar$new(total = nrow(param_grid))
# errors_case1 <- map(
#   1:nrow(param_grid),
#   ~ {
#     sim_error_case1(
#       param_grid[.x, 2],
#       param_grid[.x, 3],
#       param_grid[.x, 1],
#       param_grid[.x, 4],
#       print_plots = TRUE
#     )
#     pb$tick()
#   }
# )
errors_case1 <- future_map(
  1:nrow(param_grid),
  ~ sim_error_case1(
    param_grid[.x, 2],
    param_grid[.x, 3],
    param_grid[.x, 1],
    param_grid[.x, 4],
    print_plots = FALSE
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

errors_case2 <- future_map(
  1:nrow(param_grid),
  ~ sim_error_case2(
    param_grid[.x, 2],
    param_grid[.x, 3],
    param_grid[.x, 1],
    param_grid[.x, 4],
    print_plots = FALSE
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

errors_case3 <- future_map(
  1:nrow(param_grid),
  ~ sim_error_case3(
    param_grid[.x, 2],
    param_grid[.x, 3],
    param_grid[.x, 1],
    param_grid[.x, 4],
    print_plots = FALSE
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

errors_case4 <- future_map(
  1:nrow(param_grid),
  ~ sim_error_case4(
    param_grid[.x, 2],
    param_grid[.x, 3],
    param_grid[.x, 1],
    param_grid[.x, 4],
    print_plots = FALSE
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

errors_case5 <- future_map(
  1:nrow(param_grid),
  ~ sim_error_case5(
    param_grid[.x, 2],
    param_grid[.x, 3],
    param_grid[.x, 1],
    param_grid[.x, 4],
    print_plots = FALSE
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

errors_case6 <- future_map(
  1:nrow(param_grid),
  ~ sim_error_case6(
    param_grid[.x, 2],
    param_grid[.x, 3],
    param_grid[.x, 1],
    param_grid[.x, 4],
    print_plots = FALSE
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

errors_case7 <- future_map(
  1:nrow(param_grid),
  ~ sim_error_case7(
    param_grid[.x, 2],
    param_grid[.x, 3],
    param_grid[.x, 1],
    param_grid[.x, 4],
    print_plots = FALSE
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

errors_case8 <- future_map(
  1:nrow(param_grid),
  ~ sim_error_case8(
    param_grid[.x, 2],
    param_grid[.x, 3],
    param_grid[.x, 1],
    param_grid[.x, 4],
    print_plots = FALSE
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

errors_case9 <- future_map(
  1:nrow(param_grid),
  ~ sim_error_case9(
    param_grid[.x, 2],
    param_grid[.x, 3],
    param_grid[.x, 1],
    param_grid[.x, 4],
    print_plots = FALSE
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

errors_case10 <- future_map(
  1:nrow(param_grid),
  ~ sim_error_case10(
    param_grid[.x, 2],
    param_grid[.x, 3],
    param_grid[.x, 1],
    param_grid[.x, 4],
    print_plots = FALSE
  ),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)