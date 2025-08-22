library(doFuture)
library(progressr)
library(stringr)
library(tidyverse)

handlers(global = TRUE)
handlers("progress")

setwd("code")
source("functions/estimate_volume.R")
source("functions/simulations/estimate_volume_case8.R")


dir_case8 <- "../output/simulations/case8"
files_case8 <- list.files(dir_case8)

files_case9 <- list.files("../output/simulations/case9")

ncores <- parallel::detectCores()
plan(multisession, workers = ncores)

calc_volumes <- function(files_case8) {
  p <- progressor(along = files_case8)
  volume_dfs <- foreach(
    file_idx = seq_along(files_case8),
    .options.future = list(seed = TRUE)
  ) %dofuture%
    {
      file_name <- files_case8[file_idx]
      file_stripped <- str_replace(file_name, ".RDS", "")
      file_parsed <- str_split(file_stripped, "_")

      duration_val <- as.numeric(file_parsed[[1]][2])
      interval_val <- as.numeric(file_parsed[[1]][4]) / 100
      obs_noise_val <- as.numeric(file_parsed[[1]][7]) / 100
      amplitude_noise_val <- as.numeric(file_parsed[[1]][10]) / 100
      period_noise_val <- ifelse(
        as.numeric(file_parsed[[1]][[13]]) > 1,
        as.numeric(file_parsed[[1]][[13]]) / 100,
        as.numeric(file_parsed[[1]][[13]])
      )
      time_trend_val <- file_parsed[[1]][16]
      time_change_val <- as.numeric(file_parsed[[1]][19]) / 100

      sim <- readRDS(file.path(dir_case8, file_name))
      sim_volume <- estimate_volume_case8(
        sim,
        time_trend_val,
        time_change_val,
        n_points = 10000
      )

      n_obs <- length(sim_volume$times)
      sim_run <- rep(file_idx, n_obs)
      duration <- rep(duration_val, n_obs)
      interval <- rep(interval_val, n_obs)
      obs_noise <- rep(obs_noise_val, n_obs)
      amplitude_noise <- rep(amplitude_noise_val, n_obs)
      period_noise <- rep(period_noise_val, n_obs)
      time_trend <- rep(time_trend_val, n_obs)
      time_change <- rep(time_change_val, n_obs)

      time <- sim_volume$times
      observed_volume <- sim_volume$observed
      true_volume <- sim_volume$true
      lpme_aug_volume <- sim_volume$lpme_augmented
      pme_aug_volume <- sim_volume$pme_augmented
      lpme_partitioned_volume <- sim_volume$lpme_part
      pme_partitioned_volume <- sim_volume$pme_part
      pc_partitioned_volume <- sim_volume$pc_part
      lpme_partitioned_mesh_volume <- sim_volume$lpme_part_mesh
      pme_partitioned_mesh_volume <- sim_volume$pme_part_mesh

      vol_df <- tibble(
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
        pme_partitioned_volume = pme_partitioned_volume,
        pc_partitioned_volume = pc_partitioned_volume,
        lpme_partitioned_mesh_volume = lpme_partitioned_mesh_volume,
        pme_partitioned_mesh_volume = pme_partitioned_mesh_volume
      )
      p(sprintf("run: %g", file_idx))

      vol_df
    }
  df_full <- reduce(volume_dfs, bind_rows)
}

print("Calculating case 8 volumes")
sim_volume_df <- calc_volumes(files_case8)
write_csv(sim_volume_df, "../output/simulation_volume_case8.csv")

print("Calculating case 9 volumes")
sim_volume_df9 <- calc_volumes(files_case9)
write_csv(sim_volume_df9, "../output/simulation_volume_case9.csv")
