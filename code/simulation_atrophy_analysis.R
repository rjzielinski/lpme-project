library(tidyverse)

source("code/functions/estimate_volume.R")
source("code/functions/simulations/estimate_volume_case8.R")


dir_case8 <- "output/simulations/case8"
files_case8 <- list.files(dir_case8)

sim_run <- vector(mode = "numeric", length = 0)
duration <- vector(mode = "numeric", length = 0)
interval <- vector(mode = "numeric", length = 0)
obs_noise <- vector(mode = "numeric", length = 0)
amplitude_noise <- vector(mode = "numeric", length = 0)
period_noise <- vector(mode = "numeric", length = 0)
time_trend <- vector(mode = "character", length = 0)
time_change <- vector(mode = "numeric", length = 0)

time <- vector(mode = "numeric", length = 0)
observed_volume <- vector(mode = "numeric", length = 0)
true_volume <- vector(mode = "numeric", length = 0)
lpme_aug_volume <- vector(mode = "numeric", length = 0)
pme_aug_volume <- vector(mode = "numeric", length = 0)
lpme_partitioned_volume <- vector(mode = "numeric", length = 0)
pme_partitioned_volume <- vector(mode = "numeric", length = 0)

for (file_idx in seq_along(files_case8)) {
  file_name <- files_case8[file_idx]
  file_stripped <- str_replace(file_name, ".RDS", "")
  file_parsed <- str_split(file_stripped, "_")

  duration_val <- as.numeric(file_parsed[[1]][2])
  interval_val <- as.numeric(file_parsed[[1]][4]) / 100
  obs_noise_val <- as.numeric(file_parsed[[1]][7]) / 100
  amplitude_noise_val <- as.numeric(file_parsed[[1]][10]) / 100
  period_noise_val <- as.numeric(file_parsed[[1]][[13]])
  time_trend_val <- file_parsed[[1]][16]
  time_change_val <- as.numeric(file_parsed[[1]][19]) / 100

  sim <- readRDS(file.path(dir_case8, file_name))
  sim_volume <- estimate_volume_case8(sim, n_points = 10000)

  n_obs <- length(sim_volume$times)
  sim_run <- c(sim_run, rep(file_idx, n_obs))
  duration <- c(duration, rep(duration_val, n_obs))
  interval <- c(interval, rep(interval_val, n_obs))
  obs_noise <- c(obs_noise, rep(obs_noise_val, n_obs))
  amplitude_noise <- c(amplitude_noise, rep(amplitude_noise_val, n_obs))
  period_noise <- c(period_noise, rep(period_noise_val, n_obs))
  time_trend <- c(time_trend, rep(time_trend_val, n_obs))
  time_change <- c(time_change, rep(time_change_val, n_obs))

  time <- c(time, sim_volume$times)
  observed_volume <- c(observed_volume, sim_volume$observed)
  true_volume <- c(true_volume, sim_volume$true)
  lpme_aug_volume <- c(lpme_aug_volume, sim_volume$lpme_augmented)
  pme_aug_volume <- c(pme_aug_volume, sim_volume$pme_augmented)
  lpme_partitioned_volume <- c(lpme_partitioned_volume, sim_volume$lpme_part)
  pme_partitioned_volume <- c(pme_partitioned_volume, sim_volume$pme_part)

  cat("\rFinished", file_idx, "of", length(files_case8))
}

sim_volume_df <- tibble(
  run = sim_run,
  duration = duration,
  interval = interval,
  obs_noise = obs_noise,
  amplitude_noise = amplitude_noise,
  period_noise = period_noise,
  time_trend = time_trend,
  time_change = time_change,
  time = time,
  observed_volume = observed_volume,
  true_volume = true_volume,
  lpme_aug_volume = lpme_aug_volume,
  pme_aug_volume = pme_aug_volume,
  lpme_partitioned_volume = lpme_partitioned_volume,
  pme_partitioned_volume = pme_partitioned_volume
)

write_csv(sim_volume_df, "../output/simulation_volume_case8.csv")

