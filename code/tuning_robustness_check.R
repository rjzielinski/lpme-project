library(tidyverse)
library(plotly)
library(pracma)
library(foreach)
library(doFuture)
library(furrr)
library(here)
library(progressr)
library(pme)
library(princurve)
library(doSNOW)
library(foreach)
library(doRNG)

setwd(here("code"))

source("prinSurf_v3.R")
source("functions/simulations/simulate_data.R")
source("functions/simulations/preprocess_data.R")
source("functions/simulations/fit_models.R")
source("functions/simulations/evaluate_models.R")
source("functions/simulations/display_results.R")
source("functions/calculate_lpme_reconstructions.R")
source("functions/calculate_pme_reconstructions.R")

cores <- availableCores() - 4
plan(multicore, workers = cores)

n_sims <- 1000

set.seed(500)

with_progress({
  p <- progressor(n_sims)

  sim_results <- foreach(
    sim_idx = seq_len(n_sims),
    .options.future = list(seed = TRUE)
  ) %dofuture%
    {
      sim_data <- simulate_data(
        duration = 5,
        interval = 0.5,
        case = 1,
        obs_noise = 0.1,
        amplitude_noise = 0.1,
        period_noise = 0.1,
        time_trend = "linear",
        time_change = 0.2,
        N = 1000
      )
      sim_df <- sim_data$df

      sim_processed <- preprocess_data(sim_df, case = 1, d = 1, D = 2)

      default_gamma <- -25:5
      deltas <- rnorm(length(default_gamma), mean = 0, sd = 0.05)
      gamma_vec <- c(
        default_gamma,
        default_gamma + (abs(default_gamma) * deltas)
      )
      gamma_vec <- exp(gamma_vec)

      gamma_idx <- rep(seq_along(default_gamma), 2)

      gamma_order <- sort(gamma_vec, index.return = TRUE)
      gamma_vec <- gamma_order$x
      gamma_idx <- gamma_idx[gamma_order$ix]

      sim_lpme <- lpme(
        as.matrix(sim_processed$df_observed),
        d = 1,
        gamma = gamma_vec,
        verbose = FALSE,
        print_plots = FALSE
      )

      opt_idx <- which.min(sim_lpme$msd)
      opt_msd <- sim_lpme$msd[opt_idx]
      opt_gamma <- gamma_vec[opt_idx]
      delta_idx <- gamma_idx[opt_idx]
      delta_val <- deltas[delta_idx]

      gamma_msds <- sim_lpme$msd[gamma_idx == delta_idx]
      mod_msd <- gamma_msds[gamma_msds != opt_msd]

      msd_diff <- mod_msd - opt_msd

      p()

      length_msd <- length(sim_lpme$msd)

      sim_mat <- cbind(
        rep(sim_idx, length_msd),
        gamma_vec[1:length_msd],
        sim_lpme$msd,
        rep(opt_gamma, length_msd),
        rep(delta_val, length_msd),
        rep(opt_msd, length_msd),
        rep(mod_msd, length_msd)
      )
      sim_mat
    }
})

sim_results_df <- do.call(rbind, sim_results) %>%
  as.data.frame() %>%
  setNames(c(
    "sim_idx",
    "gamma",
    "msd",
    "opt_gamma",
    "opt_delta",
    "opt_msd",
    "mod_msd"
  ))


write_csv(sim_results_df, here("output/lpme_tuning_robustness_check.csv"))
