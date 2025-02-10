library(tidyverse)
library(plotly)
library(pracma)
library(foreach)
library(doParallel)
library(furrr)
library(progress)
library(princurve)
library(devtools)
library(plot3D)
library(latex2exp)

load_all()

source("~/Documents/brown/research/lpme-project/code/functions/sim_data.R")
source("~/Documents/brown/research/lpme-project/code/functions/calc_pme_est.R")
source("~/Documents/brown/research/lpme-project/code/functions/calc_lpme_est.R")
source("~/Documents/brown/research/lpme-project/code/prinSurf_v3.R")

lambda <- exp(-15:5)

max_time <- 1
interval <- 0.25
amp_noise <- 0.5
shape_noise <- 0.25
n <- 1000
time_change <- 0
time_trend <- "constant"

set.seed(1000)

time_vals <- seq(0, max_time, interval)
sim_list <- lapply(
  time_vals,
  sim_data,
  case = 1,
  noise = 0.15,
  amp_noise = amp_noise,
  period_noise = shape_noise,
  N = n,
  time_change = time_change,
  time_trend = time_trend
)

sim_df <- matrix(ncol = ncol(sim_list[[1]][[1]]))
true_vals <- matrix(ncol = ncol(sim_list[[1]][[2]]))
for (i in 1:length(sim_list)) {
  sim_df <- rbind(sim_df, sim_list[[i]][[1]])
  true_vals <- rbind(true_vals, sim_list[[i]][[2]])
}
sim_df <- sim_df[-1, ]
true_vals <- true_vals[-1, ]

png(
  "~/Documents/brown/research/lpme-project/paper/figures/sim_case1_plot2.png",
  width = 2500,
  height = 2500,
  res = 500
)
scatter3D(
  x = sim_df[, 1],
  y = sim_df[, 2],
  z = sim_df[, 3],
  pch = 20,
  col = alpha("black", 0.15),
  colkey = FALSE,
  # xlab = "",
  # ylab = "",
  # zlab = ""
  xlab = "Time",
  ylab = "",
  zlab = "",
  cex = 0.5
)
text3D(
  1.17, 
  -2.5, 
  0, 
  labels = expression(X[1]), 
  add = TRUE
)
text3D(
  -1.25, 
  -1.25, 
  -3, 
  labels = expression(X[2]), 
  add = TRUE
)

dev.off()

data <- sim_df
d <- 1
tuning_para_seq <- NULL
alpha <- 0.05
max_clusters <- 500
epsilon <- 0.05
max_iter <- 100
verbose <- TRUE
print_plots <- TRUE
increase_threshold <- 1.05
init <- "full"
smoothing_method <- "spline"

time_points <- unique(data[, 1])
if (is.null(tuning_para_seq)) {
  if (smoothing_method == "spline") {
    tuning_para_seq <- c(0, exp(-15:10))
  } else if (smoothing_method == "gp") {
    tuning_para_seq <- 0:10 + 0.5
  }
}

initialization <- initialize_lpme(data, init, time_points, d, alpha, max_clusters, min_comp = 10 * d)

png(
  "~/Documents/brown/research/lpme-project/paper/figures/lpme_data_reduction.png",
  width = 2500,
  height = 2500,
  res = 500
)
scatter3D(
  x = initialization$timevals,
  y = initialization$centers[, 1],
  z = initialization$centers[, 2],
  pch = 18,
  col = "red",
  colkey = FALSE,
  xlab = "Time",
  ylab = "",
  zlab = ""
)
scatter3D(
  x = sim_df[, 1],
  y = sim_df[, 2],
  z = sim_df[, 3],
  pch = 20,
  col = alpha("black", 0.15),
  colkey = FALSE,
  cex = 0.5,
  add = TRUE
  
)
text3D(
  1.17, -2.5, 0, labels = expression(X[1]), cex = 0.5, add = TRUE
)
text3D(-1.25, -1.25, -3, labels = expression(X[2]), cex = 0.5, add = TRUE)

dev.off()

init_pme_list <- fit_init_pmes(data, time_points, init, initialization, d, lambda = lambda)
splines <- merge_spline_coefs(init_pme_list, d, time_points)

spline_coefficients <- splines$coef_full
x_merged <- splines$x_test
params <- splines$params
times <- splines$times
n_knots <- splines$n_knots
lambda <- splines$lambda

init_params <- seq(min(params), max(params), length.out = 1000)
init_times <- rep(1:length(time_points), each = length(init_params))
init_param_vals <- rep(init_params, length(time_points))
init_out <- map(
  1:length(init_times),
  ~ init_pme_list$pme_results[[init_times[.x]]]$embedding_map(init_param_vals[.x])
) %>%
  reduce(rbind)

png(
  "~/Documents/brown/research/lpme-project/paper/figures/lpme_initialization.png",
  width = 2500,
  height = 2500,
  res = 500 
)
scatter3D(
  x = initialization$timevals,
  y = initialization$centers[, 1],
  z = initialization$centers[, 2],
  pch = 18,
  cex = 1.5,
  col = "red",
  xlab = "Time",
  ylab = "",
  zlab = ""
)
text3D(
  1.17, -2.5, 0, labels = expression(X[1]), add = TRUE
)
text3D(-1.25, -1.25, -3, labels = expression(X[2]), add = TRUE)
scatter3D(
  x = time_points[init_times],
  y = init_out[, 1],
  z = init_out[, 2],
  colvar = init_param_vals,
  colkey = FALSE,
  pch = 20,
  add = TRUE
)
dev.off()

D_coef <- dim(spline_coefficients)[2]
D_out <- dim(x_merged)[2]
d <- dim(params)[2]
n <- dim(x_merged)[1]
I_new <- n
t_initial <- params %>%
  as.matrix()

MSE_seq_new <- vector()
SOL_coef <- list()
TNEW_new <- list()
coefs <- list()
x_funs <- list()
functions <- list()
func_coef <- list()

inv_errors <- 1 / init_pme_list$errors
weights <- inv_errors / sum(inv_errors)

for (tuning_ind in 1:length(tuning_para_seq)) {
  if (smoothing_method == "spline") {
    f_coef_list <- compute_f_coef(
      tuning_para_seq[tuning_ind],
      diag(weights),
      t_initial,
      times,
      spline_coefficients,
      1,
      D_coef,
      4 - d,
      3
    )

    f_new <- function(t) {
      coefs <- f_coef_list$f(t[1])
      coef_mat <- matrix(coefs, n_knots + d + 1, byrow = TRUE)
      return_vec <- t(coef_mat[1:n_knots, ]) %*% etaFunc(t[-1], t_initial, 4 - d) +
        t(coef_mat[(n_knots + 1):(n_knots + d + 1), ]) %*% matrix(c(1, t[-1]), ncol = 1)
      c(t[1], return_vec)
    }
  } else if (smoothing_method == "gp") {
    # invisible(
    #   capture.output({
    #     gp <- GPFDA::gpr(
    #       response = spline_coefficients,
    #       input = times,
    #       Cov = "matern",
    #       meanModel = "t",
    #       nu = tuning_para_seq[tuning_ind]
    #     )
    #   })
    # )
    #
    # f_new <- function(t) {
    #   coefs <- GPFDA::gprPredict(
    #     train = gp,
    #     inputNew = t[1],
    #     noiseFreePred = TRUE
    #   )$pred.mean %>%
    #     as.vector()
    #   coef_mat <- matrix(coefs, n_knots + d + 1, byrow = TRUE)
    #   return_vec <- t(coef_mat[1:n_knots, ]) %*% etaFunc(t[-1], t_initial, 4 - d) +
    #     t(coef_mat[(n_knots + 1):(n_knots + d + 1), ]) %*% matrix(c(1, t[-1]), ncol = 1)
    #   c(t[1], return_vec)
    # }
  }

  updated_param <- update_parameterization(
   time_points,
    t_initial,
    x_merged,
    f_new,
    n_knots,
    d,
    d,
    4 - d
  )

  if (print_plots == TRUE) {
    plot_lpme(data, f_new, t_initial, d, D_out - 1, time_points)
  }

  # Cross-validation section
  nearest_x <- calc_nearest_x(data, x_merged)
  init_param <- calc_init_param(data, updated_param$parameterization, nearest_x)

  cv_mse <- calc_mse_cv(
    leave_one_out = TRUE,
    f = f_new,
    df = data,
    init_param = init_param,
    time_points = time_points,
    r = updated_param$parameterization[, -1],
    r_initial = t_initial,
    n_knots = n_knots,
    d = d,
    d_new2 = 1,
    D_out = D_out - 1,
    D_coef = D_coef,
    lambda = lambda,
    gamma = 4 - d,
    gamma2 = 3,
    r_full2 = times,
    w = tuning_para_seq[tuning_ind],
    smoothing_method = smoothing_method
  )

  data_n <- sapply(time_points, function(x) nrow(data[data[, 1] == x, ]))
  MSE_new <- stats::weighted.mean(cv_mse, data_n)
  MSE_seq_new[tuning_ind] <- MSE_new

  if (verbose == TRUE) {
    print_mse(tuning_para_seq[tuning_ind], MSE_new)
  }

  SOL_coef[[tuning_ind]] <- ifelse(
    smoothing_method == "spline",
    f_coef_list$sol,
    NA
  )
  TNEW_new[[tuning_ind]] <- updated_param$parameterization
  coefs[[tuning_ind]] <- params
  x_funs[[tuning_ind]] <- updated_param$embedding
  functions[[tuning_ind]] <- f_new
  func_coef[[tuning_ind]] <- ifelse(
    smoothing_method == "spline",
    f_coef_list$f,
    NA
  )

  if (tuning_ind >= 8) {
    if (!is.unsorted(MSE_seq_new[(tuning_ind - 5):tuning_ind])) {
      break
    }
  }
}

optimal_ind <- min(which(MSE_seq_new == min(MSE_seq_new)))
sol_opt <- SOL_coef[[optimal_ind]]
t_new_opt <- TNEW_new[[optimal_ind]]
coefs_opt <- coefs[[optimal_ind]]
f.optimal <- functions[[optimal_ind]]

function_list <- map(func_coef, ~ function(t) {
  coefs <- .x(t[1])
  coef_mat <- matrix(coefs, n_knots + d + 1, byrow = TRUE)
  return_vec <- t(coef_mat[1:n_knots, ]) %*% etaFunc(t[-1], t_initial, 4 - d) +
        t(coef_mat[(n_knots + 1):(n_knots + d + 1), ]) %*% matrix(c(1, t[-1]), ncol = 1)
  c(t[1], return_vec)
})

test_out <- list()
test_times <- rep(seq(min(time_vals), max(time_vals), by = 0.01), each = length(init_params))
test_params <- rep(init_params, length(seq(min(time_vals), max(time_vals), by = 0.01)))
for (i in 1:length(tuning_para_seq)) {
  test_out[[i]] <- map(
    1:length(test_times),
    ~ function_list[[i]](cbind(test_times[.x], test_params[.x])
    )
  ) %>%
    reduce(rbind)
}

png(
  "~/Documents/brown/research/lpme-project/paper/figures/lpme_fitting.png",
  width = 2500,
  height = 2500,
  res = 500
)
scatter3D(
  x = initialization$timevals,
  y = initialization$centers[, 1],
  z = initialization$centers[, 2],
  pch = 18,
  col = "red",
  cex = 1.5,
  xlab = "Time",
  ylab = "",
  zlab = ""
)
text3D(
  1.17, -2.5, 0, labels = expression(X[1]), add = TRUE
)
text3D(-1.25, -1.25, -3, labels = expression(X[2]), add = TRUE)
scatter3D(
  x = test_out[[optimal_ind]][, 1],
  y = test_out[[optimal_ind]][, 2],
  z = test_out[[optimal_ind]][, 3],
  pch = 20,
  col = alpha("navy", 0.01),
  colkey = FALSE,
  add = TRUE
)
# for (i in 1:length(tuning_para_seq)) {
#   scatter3D(
#     x = test_out[[i]][, 1],
#     y = test_out[[i]][, 2],
#     z = test_out[[i]][, 3],
#     pch = 20,
#     col = alpha("navy", 0.01),
#     colkey = FALSE,
#     add = TRUE
#   )
# }
dev.off()
