library(tidyverse)
library(plotly)
library(pracma)
library(furrr)
library(progress)
library(pme)
library(princurve)
library(doSNOW)
library(foreach)
library(doRNG)

setwd("code")

source("prinSurf_v3.R")
source("functions/simulations/simulate_data.R")
source("functions/simulations/preprocess_data.R")
source("functions/simulations/fit_models.R")
source("functions/simulations/evaluate_models.R")
source("functions/simulations/display_results.R")
source("functions/calculate_lpme_reconstructions.R")
source("functions/calculate_pme_reconstructions.R")

durations <- c(1, 2, 5)
intervals <- c(0.25, 0.5)
cases <- 8:1
obs_noises <- c(0.05, 0.1, 0.25, 0.5)
amplitude_noises <- c(0.05, 0.1, 0.25, 0.5)
period_noises <- c(0.05, 0.1, 0.25, 0.5)
time_trends <- c("constant", "linear", "quadratic")
time_changes <- c(0, 0.05, 0.1, 0.25, 0.5)


parameter_df <- expand_grid(
  duration = durations,
  interval = intervals,
  case = cases,
  obs_noise = obs_noises,
  amplitude_noise = amplitude_noises,
  period_noise = period_noises,
  time_trend = time_trends,
  time_change = time_changes
)

error_df <- tibble(
  lpme_error = vector(mode = "numeric", length = nrow(parameter_df)),
  pme_error = vector(mode = "numeric", length = nrow(parameter_df)),
  pc_error = vector(mode = "numeric", length = nrow(parameter_df)),
  data_error = vector(mode = "numeric", length = nrow(parameter_df))
)

seed_states <- list()

print("Starting Simulations...")

cl <- makeCluster(parallel::detectCores() - 1)
registerDoSNOW(cl)

registerDoRNG(611492)

pb  <- txtProgressBar(max = nrow(parameter_df), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
error_df <- foreach(
  row_idx = seq_len(nrow(parameter_df)),
  .export = c(
    "simulate_data",
    "preprocess_data",
    "fit_models",
    "evaluate_models",
    "display_results"
  ),
  .packages = c("tidyverse", "pme"),
  .combine = rbind,
  .options.snow = opts
) %dopar% {
  seed_val <- .Random.seed
  seed_states[[row_idx]] <- seed_val

  duration_value <- parameter_df$duration[row_idx]
  interval_value <- parameter_df$interval[row_idx]
  case_value <- parameter_df$case[row_idx]
  obs_noise_value <- parameter_df$obs_noise[row_idx]
  amplitude_noise_value <- parameter_df$amplitude_noise[row_idx]
  period_noise_value <- parameter_df$period_noise[row_idx]
  trend_value <- parameter_df$time_trend[row_idx]
  change_value <- parameter_df$time_change[row_idx]

  # additional stability is needed for Swiss roll and spherical cases
  # to maintain structural integrity
  if (case_value == 7) {
    amplitude_noise_value <- amplitude_noise_value / 10
    period_noise_value <- period_noise_value / 10
  } else if (case_value == 8) {
    period_noise_value <- period_noise_value / 10
  }

  D_value <- case_when(
    case_value == 1 ~ 2,
    case_value == 2 ~ 2,
    case_value == 3 ~ 2,
    case_value == 4 ~ 3,
    case_value == 5 ~ 3,
    case_value == 6 ~ 3,
    case_value == 7 ~ 3,
    case_value == 8 ~ 3
  )
  d_value <- case_when(
    case_value == 1 ~ 1,
    case_value == 2 ~ 1,
    case_value == 3 ~ 1,
    case_value == 4 ~ 1,
    case_value == 5 ~ 1,
    case_value == 6 ~ 2,
    case_value == 7 ~ 2,
    case_value == 8 ~ 2
  )

  sim_data <- simulate_data(
    duration = duration_value,
    interval = interval_value,
    case = case_value,
    obs_noise = obs_noise_value,
    amplitude_noise = amplitude_noise_value,
    period_noise = period_noise_value,
    time_trend = trend_value,
    time_change = change_value,
    N = 1000
  )
  sim_df <- sim_data$df


  sim_processed <- preprocess_data(
    sim_df,
    case = case_value,
    d = d_value,
    D = D_value
  )


  sim_models <- fit_models(
    sim_processed,
    case = case_value,
    d = d_value,
    D = D_value
  )


  model_errors <- evaluate_models(
    sim_processed,
    sim_models,
    case = case_value,
    d = d_value,
    D = D_value
  )


  result_plot <- display_results(
    sim_processed,
    sim_models,
    case = case_value,
    d = d_value,
    D = D_value
  )


  simulation_output <- list(
    data = sim_data,
    processed_data = sim_processed,
    models = sim_models,
    errors = model_errors,
    plot = result_plot,
    seed = seed_val
  )

  result_dir <- paste0("../output/simulations/case", case_value)
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }

  simulation_file <- paste0(
    "duration_",
    duration_value,
    "_interval_",
    str_pad(as.character(interval_value * 100), 3, pad = "0"),
    "_obs_noise_",
    str_pad(as.character(obs_noise_value * 100), 3, pad = "0"),
    "_amplitude_noise_",
    str_pad(as.character(amplitude_noise_value * 100), 3, pad = "0"),
    "_period_noise_",
    str_pad(as.character(period_noise_value * 100), 3, pad = "0"),
    "_time_trend_",
    trend_value,
    "_time_change_",
    str_pad(as.character(change_value * 100), 3, pad = "0"),
    ".RDS"
  )


  saveRDS(simulation_output, file.path(result_dir, simulation_file))

  unlist(model_errors)
}

stopCluster(cl)

error_df <- error_df %>%
  as_tibble()
print(error_df)

print("Simulations Complete")

simulation_results <- bind_cols(parameter_df, error_df)
write_csv(simulation_results, "../output/simulation_results.csv")
