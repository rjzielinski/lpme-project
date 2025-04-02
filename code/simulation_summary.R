library(foreach)
library(stringr)
library(tidyverse)
library(lubridate)

paths <- list.files("output/simulations", recursive = TRUE, full.names = TRUE)

sim_df <- foreach(
  idx = seq_along(paths), 
  .combine = rbind, 
  .errorhandling="remove"
) %do% {
    path_split <- str_split(paths[idx], "/")[[1]]
    case_val <- gsub(pattern = "case", replacement = "", path_split[3]) |>
      as.numeric()

    file_split <- gsub(".RDS", "", path_split[4]) |>
      str_split("_")
    file_split <- file_split[[1]]

    duration_val <- as.numeric(file_split[2])
    interval_val <- as.numeric(file_split[4]) / 100
    obs_noise_val <- as.numeric(file_split[7]) / 100
    amplitude_noise_val <- as.numeric(file_split[10]) / 100
    period_noise_val <- as.numeric(file_split[13]) / 100
    time_trend_val <- file_split[16]
    time_change_val <- as.numeric(file_split[19]) / 100
  sim_object <- readRDS(paths[idx])
  unlist(sim_object$errors)
  sim_info <- list(
      duration = duration_val,
      interval = interval_val,
      case = case_val,
      obs_noise = obs_noise_val,
      amplitude_noise = amplitude_noise_val,
      period_noise = period_noise_val,
      time_trend = time_trend_val,
      time_change = time_change_val
    )

  c(unlist(sim_info), unlist(sim_object$errors))
}

sim_df <- as_tibble(sim_df)
sim_df <- sim_df |>
  rename(
    lpme_time = lpme_time.elapsed
  ) |>
  mutate(
    duration = as.numeric(duration),
    interval = as.numeric(interval),
    case = as.factor(case),
    obs_noise = as.numeric(obs_noise),
    amplitude_noise = as.numeric(amplitude_noise),
    period_noise = as.numeric(period_noise),
    time_change = as.numeric(time_change),
    lpme_error = as.numeric(lpme_error),
    pme_error = as.numeric(pme_error),
    principal_curve_error = as.numeric(principal_curve_error),
    lpme_part_error = as.numeric(lpme_part_error),
    pme_part_error = as.numeric(pme_part_error),
    principal_curve_part_error = as.numeric(principal_curve_part_error),
    data_error = as.numeric(data_error),
    lpme_time = as.numeric(lpme_time),
    pme_time = as.numeric(pme_time),
    principal_curve_time = as.numeric(principal_curve_time),
    lpme_part_time = as.numeric(lpme_part_time),
    pme_part_time = as.numeric(pme_part_time),
    principal_curve_part_time = as.numeric(principal_curve_part_time)
  )

write_csv(sim_df, "output/simulation_results.csv")
