









#| label: chunk_setup
knitr::opts_chunk$set(message = FALSE, warning = FALSE)



#| label: load_packages
library(tidyverse)
library(plotly)
library(pracma)
library(foreach)
library(doParallel)
library(furrr)
library(progress)
library(pme)
library(princurve)
library(Rfast)

source("prinSurf_v3.R")
source("functions/simulations/simulate_data.R")
source("functions/simulations/preprocess_data.R")
source("functions/simulations/fit_models.R")
source("functions/simulations/evaluate_models.R")
source("functions/simulations/display_results.R")
source("functions/calculate_lpme_reconstructions.R")
source("functions/calculate_pme_reconstructions.R")
source("functions/estimate_volume.R")
source("functions/interior_identification.R")








#| label: plot_case1
set.seed(100)
sim_data <- simulate_data(
  duration = 1,
  interval = 0.25,
  case = 1,
  obs_noise = 0.1,
  amplitude_noise = 0.25,
  period_noise = 0.05,
  time_trend = "constant",
  time_change = 0,
  N = 1000
)
sim_df <- sim_data$df

plot_ly(
  data = sim_df,
  x = ~X1,
  y = ~time,
  z = ~X2,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = "black"),
  name = "Data"
) |>
  add_markers(
    x = ~X_true1,
    y = ~time,
    z = ~X_true2,
    marker = list(size = 3, color = "red"),
    name = "True Manifold"
  )




#| label: plot_case2
set.seed(200)
sim_data <- simulate_data(
  duration = 1,
  interval = 0.25,
  case = 2,
  obs_noise = 0.1,
  amplitude_noise = 0.25,
  period_noise = 0.05,
  time_trend = "constant",
  time_change = 0,
  N = 1000
)
sim_df <- sim_data$df

plot_ly(
  data = sim_df,
  x = ~X1,
  y = ~time,
  z = ~X2,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = "black"),
  name = "Data"
) |>
  add_markers(
    x = ~X_true1,
    y = ~time,
    z = ~X_true2,
    marker = list(size = 3, color = "red"),
    name = "True Manifold"
  )





#| label: plot_case3
set.seed(300)
sim_df <- simulate_data(
  duration = 1,
  interval = 0.25,
  case = 3,
  obs_noise = 0.15,
  amplitude_noise = 0.25,
  period_noise = 0.10,
  time_trend = "constant",
  time_change = 0,
  N = 1000
)
sim_df <- sim_data$df 

plot_ly(
  data = sim_df,
  x = ~X1,
  y = ~time,
  z = ~X2,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = "black"),
  name = "Data"
) |>
  add_markers(
    x = ~X_true1,
    y = ~time,
    z = ~X_true2,
    marker = list(size = 3, color = "red"),
    name = "True Manifold"
  )




#| label: plot_case4
set.seed(400)
sim_data <- simulate_data(
  duration = 1,
  interval = 0.25,
  case = 4,
  obs_noise = 0.15,
  amplitude_noise = 0.25,
  period_noise = 0.10,
  time_trend = "constant",
  time_change = 0,
  N = 1000
)
sim_df <- sim_data$df

plot_ly(
  data = sim_df,
  x = ~X1,
  y = ~X2,
  z = ~X3,
  frame = ~time,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = "black"),
  name = "Data"
) |>
  add_markers(
    x = ~X_true1,
    y = ~X_true2,
    z = ~X_true3,
    frame = ~time,
    marker = list(size = 3, color = "red"),
    name = "True Manifold"
  )





#| label: plot_case5
set.seed(500)
sim_data <- simulate_data(
  duration = 1,
  interval = 0.25,
  case = 5,
  obs_noise = 0.15,
  amplitude_noise = 0.25,
  period_noise = 0.10,
  time_trend = "constant",
  time_change = 0,
  N = 1000
)
sim_df <- sim_data$df

plot_ly(
  data = sim_df,
  x = ~X1,
  y = ~X2,
  z = ~X3,
  frame = ~time,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = "black"),
  name = "Data"
) |>
  add_markers(
    x = ~X_true1,
    y = ~X_true2,
    z = ~X_true3,
    frame = ~time,
    marker = list(size = 3, color = "red"),
    name = "True Manifold"
  )






#| label: plot_case6
set.seed(600)
sim_data <- simulate_data(
  duration = 1,
  interval = 0.25,
  case = 6,
  obs_noise = 0.15,
  amplitude_noise = 0.25,
  period_noise = 0.10,
  time_trend = "constant",
  time_change = 0,
  N = 1000
)
sim_df <- sim_data$df

plot_ly(
  data = sim_df,
  x = ~X1,
  y = ~X2,
  z = ~X3,
  frame = ~time,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = "black"),
  name = "Data"
) |>
  add_markers(
    x = ~X_true1,
    y = ~X_true2,
    z = ~X_true3,
    frame = ~time,
    marker = list(size = 3, color = "red"),
    name = "True Manifold"
  )






#| label: plot_case7
set.seed(700)
sim_data <- simulate_data(
  duration = 1,
  interval = 0.25,
  case = 7,
  obs_noise = 0.15,
  amplitude_noise = 0.05,
  period_noise = 0.01,
  time_trend = "constant",
  time_change = 0,
  N = 1000
)
sim_df <- sim_data$df

plot_ly(
  data = sim_df,
  x = ~X1,
  y = ~X2,
  z = ~X3,
  frame = ~time,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = "black"),
  name = "Data"
) |>
  add_markers(
    x = ~X_true1,
    y = ~X_true2,
    z = ~X_true3,
    frame = ~time,
    marker = list(size = 3, color = "red"),
    name = "True Manifold"
  )






#| label: plot_case8
set.seed(800)
sim_data <- simulate_data(
  duration = 1,
  interval = 0.25,
  case = 8,
  obs_noise = 0.15,
  amplitude_noise = 0.25,
  period_noise = 0.05,
  time_trend = "constant",
  time_change = 0,
  N = 1000
)
sim_df <- sim_data$df

plot_ly(
  data = sim_df,
  x = ~X1,
  y = ~X2,
  z = ~X3,
  frame = ~time,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = "black"),
  name = "Data"
) |>
  add_markers(
    x = ~X_true1,
    y = ~X_true2,
    z = ~X_true3,
    frame = ~time,
    marker = list(size = 3, color = "red"),
    name = "True Manifold"
  )



































#| label: demo_case1
set.seed(100)
sim_data <- simulate_data(
  duration = 1,
  interval = 0.25,
  case = 1,
  obs_noise = 0.1,
  amplitude_noise = 0.25,
  period_noise = 0.05,
  time_trend = "constant",
  time_change = 0,
  N = 1000
)
sim_df <- sim_data$df

sim_processed <- preprocess_data(sim_df, case = 1, d = 1, D = 2)
sim_models <- fit_models(sim_processed, case = 1, d = 1, D = 2)
model_errors <- evaluate_models(
  sim_processed,
  sim_models,
  case = 1,
  d = 1,
  D = 2
)
model_errors
result_plot <- display_results(
  sim_processed,
  sim_models,
  case = 1,
  d = 1,
  D = 2
)
result_plot




#| label: demo_case3
set.seed(300)
sim_data <- simulate_data(
  duration = 1,
  interval = 0.25,
  case = 3,
  obs_noise = 0.1,
  amplitude_noise = 0.25,
  period_noise = 0.05,
  time_trend = "constant",
  time_change = 0,
  N = 1000
)
sim_df <- sim_data$df

sim_processed <- preprocess_data(sim_df, case = 3, d = 1, D = 2)
sim_models <- fit_models(sim_processed, case = 3, d = 1, D = 2)
model_errors <- evaluate_models(
  sim_processed,
  sim_models,
  case = 3,
  d = 1,
  D = 2
)
model_errors

result_plot <- display_results(
  sim_processed,
  sim_models,
  case = 3,
  d = 1,
  D = 2
)
result_plot




#| label: demo_case7
set.seed(700)
sim_data <- simulate_data(
  duration = 1,
  interval = 0.25,
  case = 7,
  obs_noise = 0.1,
  amplitude_noise = 0.01,
  period_noise = 0.01,
  time_trend = "constant",
  time_change = 0,
  N = 1000
)
sim_df <- sim_data$df

sim_processed <- preprocess_data(sim_df, case = 7, d = 2, D = 3)
sim_models <- fit_models(sim_processed, case = 7, d = 2, D = 3)
model_errors <- evaluate_models(
  sim_processed,
  sim_models,
  case = 7,
  d = 2,
  D = 3
)
model_errors

result_plot <- display_results(
  sim_processed,
  sim_models,
  case = 7,
  d = 2,
  D = 3
)
result_plot




#| label: demo_case8
set.seed(800)
sim_data <- simulate_data(
  duration = 1,
  interval = 0.25,
  case = 8,
  obs_noise = 0.1,
  amplitude_noise = 0.25,
  period_noise = 0.05,
  time_trend = "constant",
  time_change = 0,
  N = 1000
)
sim_df <- sim_data$df

sim_processed <- preprocess_data(sim_df, case = 8, d = 2, D = 3)
sim_models <- fit_models(sim_processed, case = 8, d = 2, D = 3)
model_errors <- evaluate_models(
  sim_processed,
  sim_models,
  case = 8,
  d = 2,
  D = 3
)
model_errors

result_plot <- display_results(
  sim_processed,
  sim_models,
  case = 8,
  d = 2,
  D = 3
)
result_plot







#| label: simulation_parameters

# set possible values of simulation parameters as discussed above
durations <- c(1, 2, 5)
intervals <- c(0.25, 0.5)
cases <- 1:8
obs_noises <- c(0.05, 0.1, 0.25, 0.5)
amplitude_noises <- c(0.05, 0.1, 0.25, 0.5)
period_noises <- c(0.05, 0.1, 0.25, 0.5)
time_trends <- c("constant", "linear", "quadratic")
time_changes <- c(0, 0.05, 0.1, 0.25, 0.5)


# create dataframe of all combinations of parameter values
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

# create empty dataframe to store error estimates
error_df <- tibble(
  lpme_error = vector(mode = "numeric", length = nrow(parameter_df)),
  pme_error = vector(mode = "numeric", length = nrow(parameter_df)),
  pc_error = vector(mode = "numeric", length = nrow(parameter_df)),
  data_error = vector(mode = "numeric", length = nrow(parameter_df))
)

# save the seed state that was set when each simulation was run
seed_states <- list()




#| label: run_simulations
print("Starting Simulations...")

cl <- makeCluster(parallel::detectCores() - 1)
# unclear why, but registerDoSNOW worked better when running with slurm
registerDoSNOW(cl)

# the registerDoRNG function allows reproducibility
# of the simulations when run in parallel
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

  # following the same analysis pipeline as described above
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



#| label: save_simulation_results
error_df <- error_df %>%
  as_tibble()
print(error_df)

print("Simulations Complete")

simulation_results <- bind_cols(parameter_df, error_df)
write_csv(simulation_results, "../output/simulation_results.csv")





#| label: sample_result
result_case1 <- readRDS(
  file.path(
    "../output/simulations/case1",
    "duration_1_interval_025_obs_noise_005_amplitude_noise_010_period_noise_010_time_trend_linear_time_change_010.RDS"
    # "duration_1_interval_050_obs_noise_010_amplitude_noise_025_period_noise_010_time_trend_linear_time_change_010.RDS"
    # "duration_1_interval_025_obs_noise_005_amplitude_noise_005_period_noise_005_time_trend_constant_time_change_025.RDS"
  )
)

case1_observed <- result_case1$processed_data$df_observed
case1_true <- result_case1$processed_data$df_true
case1_sim <- result_case1$data$df

lpme_data <- result_case1$models$lpme$reconstructions
pme_data <- result_case1$models$pme$reconstructions
pc_data <- result_case1$models$pc$reconstructions

case1_names <- c("time", "X1", "X2")
names(case1_observed) <- case1_names
names(case1_true) <- case1_names
case1_observed$source <- "Data"
case1_observed$param <- case1_sim$r1
case1_true$source <- "True Manifold"
case1_true$param <- case1_sim$r1

lpme_df <- as_tibble(lpme_data, ..name_repair = NULL)
names(lpme_df) <- case1_names
lpme_df$source <- "LPME"
lpme_df$param <- case1_sim$r1

pme_df <- as_tibble(pme_data, ..name_repair = NULL)
names(pme_df) <- case1_names
pme_df$source <- "PME"
pme_df$param <- case1_sim$r1

pc_df <- as_tibble(pc_data, ..name_repair = NULL)
names(pc_df) <- case1_names
pc_df$source <- "PC"
pc_df$param <- case1_sim$r1

result_case1$plot





#| label: result_plot
case1_data <- bind_rows(
  case1_observed,
  case1_true,
  lpme_df,
  pme_df,
  pc_df
) |>
  arrange(param)

color_values <- c("black", "blue", "purple", "green", "red")

set.seed(8195)
# subsample data to improve clarity of visualization
data_sampled <- case1_data |>
  group_by(time, source) |>
  sample_n(100)



#| label: result_facet_plot
ggplot(case1_data, aes(x = X1, y = X2)) +
  facet_wrap(~time, nrow = 3) +
  geom_point(
    data = subset(case1_data, source == "Data"),
    color = "black",
    alpha = 0.2,
    size = 0.25
  ) +
  geom_point(
    data = subset(case1_data, source != "Data"),
    aes(color = source, shape = source),
    size = 0.2
  ) +
  scale_color_manual(values = color_values[-1])




result_case1$errors











dir_case8 <- "../output/simulations/case8"
files_case8 <- list.files(dir_case8)

file_name <- files_case8[1]
file_stripped <- str_replace(file_name, ".RDS", "")
file_parsed <- str_split(file_stripped, "_")

duration_val <- as.numeric(file_parsed[[1]][2])
interval_val <- as.numeric(file_parsed[[1]][4]) / 100
time_trend_val <- file_parsed[[1]][16]
time_change_val <- as.numeric(file_parsed[[1]][19]) / 100

sim <- readRDS(file.path(dir_case8, file_name))

sim_data <- sim$processed_data$df
time_points <- sim_data |>
  select(time) |>
  unlist() |>
  unique()

sim_mat <- sim$processed_data$df_observed |>
  select(time, contains("X")) |>
  as.matrix()

radii <- map(
  seq_along(sim$data$amplitude_values),
  ~ sim$data$amplitude_values[[.x]][1]
) |>
  reduce(c)

observed_volumes <- (4 / 3) * pi * radii^3

time_adjustments <- case_when(
  time_trend_val == "constant" ~ 0 * time_points,
  time_trend_val == "linear" ~ time_points,
  time_trend_val == "quadratic" ~ time_points^2
)

# scale time adjustments to be proportion of 1
if (max(time_adjustments) == 0) {
  time_adjustments_scaled <- time_adjustments * time_change_val * 0
} else {
  time_adjustments_scaled <- (
    time_adjustments / max(time_adjustments)
  ) * time_change_val
}
# subtract from initial amplitude to replicate effects of atrophy
true_radii <- 1 - time_adjustments_scaled
true_volumes <- (4 / 3) * pi * true_radii^3



lpme_interior_points <- matrix(ncol = ncol(sim_mat))
lpme_volume <- vector(mode = "numeric", length = length(time_points))
pme_interior_points <- matrix(ncol = ncol(sim_mat))
pme_volumes <- vector(mode = "numeric", length = length(time_points))

df_observed <- sim$processed_data$df_observed |>
  select(time, contains("X")) |>
  as.matrix()
df_true <- sim$processed_data$df_true |>
  select(time, contains("X")) |>
  as.matrix()
df_lpme <- sim$models$lpme$reconstructions[, 1:4]
df_pme <- sim$models$pme$reconstructions[, 1:4]

for (time_idx in seq_along(time_points)) {
  temp_observed <- df_observed[df_observed[, 1] == time_points[time_idx], -1]
  temp_true <- df_true[df_true[, 1] == time_points[time_idx], -1]
  temp_lpme <- df_lpme[df_lpme[, 1] == time_points[time_idx], -1]
  temp_pme <- df_pme[df_pme[, 1] == time_points[time_idx], -1]

  scaling_factor <- 0.01
  temp_observed_scaled <- round(temp_observed / scaling_factor)
  temp_true_scaled <- round(temp_true / scaling_factor)
  temp_lpme_scaled <- round(temp_lpme / scaling_factor)
  temp_pme_scaled <- round(temp_pme / scaling_factor)
}




lpme_interior_points <- matrix(ncol = ncol(sim_mat))
lpme_volumes <- vector(mode = "numeric", length = length(time_points))

pme_interior_points <- matrix(ncol = ncol(sim_mat))
pme_volumes <- vector(mode = "numeric", length = length(time_points))

for (time_idx in seq_along(time_points)) {
  temp_data <- sim_mat[sim_mat[, 1] == time_points[time_idx], -1]

  data_limits <- colMinsMaxs(temp_data)

  limit_scaler <- 0.1
  x_seq <- seq(
    data_limits[1, 1] - abs(data_limits[1, 1] * limit_scaler),
    data_limits[2, 1] + abs(data_limits[2, 1] * limit_scaler),
    by = 0.1
  )
  y_seq <- seq(
    data_limits[1, 2] - abs(data_limits[1, 2] * limit_scaler),
    data_limits[2, 2] + abs(data_limits[2, 2] * limit_scaler),
    by = 0.1
  )
  z_seq <- seq(
    data_limits[1, 3] - abs(data_limits[1, 3] * limit_scaler),
    data_limits[2, 3] + abs(data_limits[2, 3] * limit_scaler),
    by = 0.1
  )

  point_volume <- 0.1^3

  lpme_candidates  <- expand.grid(x = x_seq, y = y_seq, z = z_seq) |>
    as.matrix()
  pme_candidates <- expand.grid(x = x_seq, y = y_seq, z = z_seq) |>
    as.matrix()

  temp_pme <- sim$models$pme$pme[[time_idx]]
  temp_pme_embedding <- temp_pme$embedding_map
  temp_pme_coefs <- temp_pme$coefs[[which.min(temp_pme$MSD)]]
  temp_pme_params <- temp_pme$parameterization[[which.min(temp_pme$MSD)]]

  lpme_est <- sim$models$lpme$lpme
  temp_lpme_coefs <- lpme_est$sol_coef_functions[[
    which.min(lpme_est$msd)
  ]](time_points[time_idx])
  temp_lpme_params <- lpme_est$parameterization_list[[which.min(lpme_est$msd)]]
  temp_lpme_params <- temp_lpme_params[
    temp_lpme_params[, 1] == time_points[time_idx],
    -1
  ]
  lpme_n_knots <- nrow(temp_lpme_params)
  d <- ncol(temp_lpme_params)
  coef_mat <- matrix(
    temp_lpme_coefs,
    lpme_n_knots + d + 1,
    byrow = TRUE
  )

  temp_lpme_embedding <- function(r) {
    t(coef_mat[1:lpme_n_knots, ]) %*% pme::etaFunc(r, temp_lpme_params, 4 - d) +
      t(coef_mat[(lpme_n_knots + 1):(lpme_n_knots + d + 1), ]) %*%
        matrix(c(1, r), ncol = 1)
  }

  lpme_interior_voxels <- interior_identification(
    temp_lpme_embedding,
    coef_mat,
    temp_lpme_params,
    lpme_candidates,
    c(0, 0, 0)
  )
  pme_interior_voxels <- interior_identification(
    temp_pme_embedding,
    temp_pme_coefs,
    temp_pme_params,
    pme_candidates,
    c(0, 0, 0)
  )

  lpme_interior <- cbind(
    time_points[time_idx],
    lpme_candidates[lpme_interior_voxels, ]
  )
  lpme_interior_points <- rbind(lpme_interior_points, lpme_interior)
  lpme_volumes[time_idx] <- nrow(lpme_interior) * point_volume

  pme_interior <- cbind(
    time_points[time_idx],
    pme_candidates[pme_interior_voxels, ]
  )
  pme_interior_points <- rbind(pme_interior_points, pme_interior)
  pme_volumes[time_idx] <- nrow(pme_interior) * point_volume
}
