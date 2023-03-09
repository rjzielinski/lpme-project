library(tidyverse)
library(plotly)
library(pracma)
library(profvis)
library(foreach)
library(doParallel)

source("code/functions/sim_data.R")
source("code/functions/calc_pme_est.R")
source("code/functions/calc_lpme_est.R")
source("code/pme.R")
source("code/lpme.R")



### Simulation Case 1

noise_vals <- seq(0, 2, 0.1)
max_times <- c(1, 2, 5, 8, 10)
intervals <- c(1 / 12, 0.25, 0.5, 1)

param_grid <- expand.grid(max_times, intervals, noise_vals)

# ncores <- parallel::detectCores() / 2
# cl <- parallel::makeCluster(ncores, type = "FORK")
# doParallel::registerDoParallel(cl = cl)


set.seed(1286)
foreach (idx = 1:nrow(param_grid)) %do% {
# for (time_idx in 1:length(max_times)) {
  time_vals <- seq(0, param_grid[idx, 1], param_grid[idx, 2])
  sim_df <- lapply(
    time_vals,
    sim_data,
    case = 1,
    noise = 0.15,
    shape_noise = param_grid[idx, 3]
  ) %>%
    reduce(rbind)

  # lpme_result <- lpme(sim_df, 1, print_plots = FALSE)
  lpme_result <- lpme(sim_df, 1)
  lpme_vals <- calc_lpme_est(lpme_result, sim_df)
  pme_result <- list()
  pme_vals <- list()
  for (t in 1:length(time_vals)) {
    temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
    pme_result[[t]] <- pme(temp_data, d = 1)
    pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
  }
  pme_vals <- reduce(pme_vals, rbind)

  tau <- sim_df[, 2]
  true_vals <- matrix(nrow = nrow(sim_df), ncol = ncol(sim_df))
  true_vals[, 1] <- sim_df[, 1]
  true_vals[, 2] <- tau
  true_vals[, 3] <- sin(tau + (pi / 2))

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
    times = time_vals,
    noise = param_grid[idx, 3],
    lpme_result = lpme_result,
    pme_results = pme_results,
    lpme_error = lpme_error,
    pme_error = pme_error
  )

  saveRDS(
    sim_case1,
    paste0(
      "simulations/case1/duration_",
      str_pad(as.character(param_grid[idx, 1]), 2, side = "left", pad = "0"),
      "_interval_",
      as.character(100 * param_grid[idx, 2]),
      "_noise_",
      as.character(100 * param_grid[idx, 3]),
      ".RDS"
    )
  )
}
parallel::stopCluster(cl = cl)

### Simulation Case 2

time_vals <- seq(0, 10, 2)

set.seed(100)
df_list <- lapply(
  time_vals,
  sim_D2d1_case2,
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  noise = 0.15,
  time_noise = 1,
  shape_noise = 0.3
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 1)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 0.1)
grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) &
    (sim_pred[, dim_idx] < idx_max)
}

r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2][r_inrange, 1]))
r_max <- max(unlist(grid_mat[, 2][r_inrange, 1]))
r_vals <- seq(
  r_min,
  r_max,
  0.1
)

grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r", "x", "y")

plot_ly(
  sim_pred_full_df,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "markers"
)

### Simulation Case 3

time_vals <- seq(0, 10, 2)

set.seed(100)
df_list <- lapply(
  time_vals,
  sim_D2d1_case3,
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  noise = 0.15,
  time_noise = 0.5,
  shape_noise = 0.3
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 1)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 0.1)
grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) &
    (sim_pred[, dim_idx] < idx_max)
}

r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2][r_inrange, 1]))
r_max <- max(unlist(grid_mat[, 2][r_inrange, 1]))
r_vals <- seq(
  r_min,
  r_max,
  0.1
)

grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r", "x", "y")

plot_ly(
  sim_pred_full_df,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "markers"
)

## Simulation Case 4

time_vals <- seq(0, 1, 0.5)

set.seed(100)
df_list <- lapply(
  time_vals,
  sim_D2d1_case4,
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  noise = 0.15,
  time_noise = 0.5,
  shape_noise = 0.3
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 1)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 0.1)
grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) &
    (sim_pred[, dim_idx] < idx_max)
}

r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2][r_inrange, 1]))
r_max <- max(unlist(grid_mat[, 2][r_inrange, 1]))
r_vals <- seq(
  r_min,
  r_max,
  0.1
)

grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r", "x", "y")

plot_ly(
  sim_pred_full_df,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "markers"
)

### D3, d1 Simulation Case 1

time_vals <- seq(0, 10, 2)

set.seed(100)
df_list <- lapply(
  time_vals,
  sim_D3d1_case1,
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  depth_multiplier = 1,
  noise = 0.15,
  time_noise = 1
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 1)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 0.1)
grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) &
    (sim_pred[, dim_idx] < idx_max)
}

r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2][r_inrange, 1]))
r_max <- max(unlist(grid_mat[, 2][r_inrange, 1]))
if (sum(r_inrange) == 0) {
  r_min <- -10
  r_max <- 10
}
r_vals <- seq(
  r_min,
  r_max,
  0.1
)

grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r", "x", "y")

plot_ly(
  sim_pred_full_df,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "markers"
)

### D3, d1 Simulation Case 2

time_vals <- seq(0, 10, 2)

set.seed(100)
df_list <- lapply(
  time_vals,
  sim_D3d1_case2,
  vertical_multiplier = 0,
  horizontal_multiplier = 0,
  depth_multiplier = 0,
  noise = 0.15,
  time_noise = 1
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 1)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 0.1)
grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(data_points) - 1)
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) &
    (sim_pred[, dim_idx] < idx_max)
}

r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2][r_inrange, 1]))
r_max <- max(unlist(grid_mat[, 2][r_inrange, 1]))
if (sum(r_inrange) == 0) {
  r_min <- -10
  r_max <- 10
}
r_vals <- seq(
  r_min,
  r_max,
  0.1
)

grid_mat <- expand_grid(time_vals, r_vals)

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(data_points) - 1)
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r", "x", "y", "z")

plot_ly(
  sim_pred_full_df,
  x = ~x,
  y = ~y,
  z = ~z,
  frame = ~time,
  opacity = 0.5,
  type = "scatter3d",
  mode = "markers"
)

### D3, d2 Simulation Case 1

time_vals <- seq(0, 10, 2)

set.seed(100)
df_list <- lapply(
  time_vals,
  sim_D3d2_case1,
  vertical_multiplier = 0,
  horizontal_multiplier = 0,
  depth_multiplier = 0,
  noise = 0.15,
  time_noise = 0.25,
  shape_noise = 0.2
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 2)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 1)

r_list <- lapply(numeric(2), function(x) r_vals)
r_mat <- as.matrix(expand.grid(r_list))

grid_mat <- expand_grid(time_vals, r_mat) %>%
  as.matrix()

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) &
    (sim_pred[, dim_idx] < idx_max)
}

r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2:3][r_inrange,]))
r_max <- max(unlist(grid_mat[, 2:3][r_inrange,]))
if (sum(r_inrange) == 0) {
  r_min <- -10
  r_max <- 10
}
r_vals <- seq(
  r_min,
  r_max,
  0.2
)
r_list <- lapply(numeric(2), function(x) r_vals)
r_mat <- as.matrix(expand.grid(r_list))

grid_mat <- expand_grid(time_vals, r_mat) %>%
  as.matrix()

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
  progress(i, nrow(sim_pred))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r1", "r2", "x1", "x2", "x3")

range_x1 <- max(sim_pred_full_df$x1) - min(sim_pred_full_df$x1)
range_x2 <- max(sim_pred_full_df$x2) - min(sim_pred_full_df$x2)
range_x3 <- max(sim_pred_full_df$x3) - min(sim_pred_full_df$x3)

plot_ly(
  sim_pred_full_df,
  x = ~x1,
  y = ~x2,
  z = ~x3,
  frame = ~time,
  opacity = 0.5,
  type = "scatter3d",
  mode = "markers"
) %>%
  layout(
    scene = list(
      xaxis = list(
        range = c(
          min(sim_pred_full_df$x1) - (0.1 * range_x1),
          max(sim_pred_full_df$x1) + (0.1 * range_x1)
        )
      ),
      yaxis = list(
        range = c(
          min(sim_pred_full_df$x2) - (0.1 * range_x2),
          max(sim_pred_full_df$x2) + (0.1 * range_x2)
        )
      ),
      zaxis = list(
        range = c(
          min(sim_pred_full_df$x3) - (0.1 * range_x3),
          max(sim_pred_full_df$x3) + (0.1 * range_x3)
        )
      )
    )
  )

### D3, d2 Simulation Case 2

time_vals <- 0:5

set.seed(100)
df_list <- lapply(
  time_vals,
  sim_D3d2_case2,
  vertical_multiplier = 0.1,
  horizontal_multiplier = 0.1,
  depth_multiplier = 0.1,
  noise = 0.15,
  time_noise = 0.1,
  shape_noise = 0.0005,
  curvature = 0.5
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 2)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 1)

r_list <- lapply(numeric(2), function(x) r_vals)
r_mat <- as.matrix(expand.grid(r_list))

grid_mat <- expand_grid(time_vals, r_mat) %>%
  as.matrix()

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) &
    (sim_pred[, dim_idx] < idx_max)
}

r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2:3][r_inrange,]))
r_max <- max(unlist(grid_mat[, 2:3][r_inrange,]))
if (sum(r_inrange) == 0) {
  r_min <- -10
  r_max <- 10
}
r_vals <- seq(
  r_min,
  r_max,
  0.2
)
r_list <- lapply(numeric(2), function(x) r_vals)
r_mat <- as.matrix(expand.grid(r_list))

grid_mat <- expand_grid(time_vals, r_mat) %>%
  as.matrix()

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
  progress(i, nrow(sim_pred))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r1", "r2", "x1", "x2", "x3")

range_x1 <- max(sim_pred_full_df$x1) - min(sim_pred_full_df$x1)
range_x2 <- max(sim_pred_full_df$x2) - min(sim_pred_full_df$x2)
range_x3 <- max(sim_pred_full_df$x3) - min(sim_pred_full_df$x3)

plot_ly(
  sim_pred_full_df,
  x = ~x1,
  y = ~x2,
  z = ~x3,
  frame = ~time,
  opacity = 0.5,
  type = "scatter3d",
  mode = "markers"
) %>%
  layout(
    scene = list(
      xaxis = list(
        range = c(
          min(sim_pred_full_df$x1) - (0.1 * range_x1),
          max(sim_pred_full_df$x1) + (0.1 * range_x1)
        )
      ),
      yaxis = list(
        range = c(
          min(sim_pred_full_df$x2) - (0.1 * range_x2),
          max(sim_pred_full_df$x2) + (0.1 * range_x2)
        )
      ),
      zaxis = list(
        range = c(
          min(sim_pred_full_df$x3) - (0.1 * range_x3),
          max(sim_pred_full_df$x3) + (0.1 * range_x3)
        )
      )
    )
  )

## D4, d2, Case 1

time_vals <- 0:5

set.seed(100)
df_list <- lapply(
  time_vals,
  sim_D5d2_case2,
  vertical_multiplier = 0.1,
  horizontal_multiplier = 0.1,
  depth_multiplier = 0.1,
  noise = 0.15,
  time_noise = 0.1,
  shape_noise = 0.0001
  # curvature = 1
)

data_points <- reduce(df_list, rbind)[, -(5:6)]
data_points_sph <- cart2sph(data_points[, 2:4])
data_points <- cbind(data_points, data_points_sph)

test_df <- data_points[data_points[, 1] == 0, -1]
pme_result <- pme(test_df, 2, print_plots = TRUE)

sim_result <- lpme(data_points, 2)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 1)

r_list <- lapply(numeric(2), function(x) r_vals)
r_mat <- as.matrix(expand.grid(r_list))

grid_mat <- expand_grid(time_vals, r_mat) %>%
  as.matrix()

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) &
    (sim_pred[, dim_idx] < idx_max)
}

r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2:3][r_inrange,]))
r_max <- max(unlist(grid_mat[, 2:3][r_inrange,]))
if (sum(r_inrange) == 0) {
  r_min <- -10
  r_max <- 10
}
r_vals <- seq(
  r_min,
  r_max,
  0.2
)
r_list <- lapply(numeric(2), function(x) r_vals)
r_mat <- as.matrix(expand.grid(r_list))

grid_mat <- expand_grid(time_vals, r_mat) %>%
  as.matrix()

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
  progress(i, nrow(sim_pred))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r1", "r2", "x1", "x2", "x3")

range_x1 <- max(sim_pred_full_df$x1) - min(sim_pred_full_df$x1)
range_x2 <- max(sim_pred_full_df$x2) - min(sim_pred_full_df$x2)
range_x3 <- max(sim_pred_full_df$x3) - min(sim_pred_full_df$x3)

plot_ly(
  sim_pred_full_df,
  x = ~x1,
  y = ~x2,
  z = ~x3,
  frame = ~time,
  opacity = 0.5,
  type = "scatter3d",
  mode = "markers"
) %>%
  layout(
    scene = list(
      xaxis = list(
        range = c(
          min(sim_pred_full_df$x1) - (0.1 * range_x1),
          max(sim_pred_full_df$x1) + (0.1 * range_x1)
        )
      ),
      yaxis = list(
        range = c(
          min(sim_pred_full_df$x2) - (0.1 * range_x2),
          max(sim_pred_full_df$x2) + (0.1 * range_x2)
        )
      ),
      zaxis = list(
        range = c(
          min(sim_pred_full_df$x3) - (0.1 * range_x3),
          max(sim_pred_full_df$x3) + (0.1 * range_x3)
        )
      )
    )
  )

## D4, d2, Case 2

time_vals <- 0:5

set.seed(100)
df_list <- lapply(
  time_vals,
  sim_D4d2_case2,
  vertical_multiplier = 0.1,
  horizontal_multiplier = 0.1,
  depth_multiplier = 0.1,
  noise = 0.15,
  time_noise = 0.1,
  shape_noise = 0.0001,
  curvature = 1
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 2)

time_vals <- seq(0, 5, 0.1)
r_vals <- seq(-10, 10, 1)

r_list <- lapply(numeric(2), function(x) r_vals)
r_mat <- as.matrix(expand.grid(r_list))

grid_mat <- expand_grid(time_vals, r_mat) %>%
  as.matrix()

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
}

idx_inrange <- matrix(nrow = dim(sim_pred)[1], ncol = dim(sim_pred)[2])
for (dim_idx in 1:dim(sim_pred)[2]) {
  idx_range <- max(data_points[, dim_idx + 1]) - min(data_points[, dim_idx + 1])
  idx_min <- min(data_points[, dim_idx + 1]) - (0.2 * idx_range)
  idx_max <- max(data_points[, dim_idx + 1]) + (0.2 * idx_range)
  idx_inrange[, dim_idx] <- (sim_pred[, dim_idx] > idx_min) &
    (sim_pred[, dim_idx] < idx_max)
}

r_inrange <- rowSums(idx_inrange) == dim(sim_pred)[2]
r_min <- min(unlist(grid_mat[, 2:3][r_inrange,]))
r_max <- max(unlist(grid_mat[, 2:3][r_inrange,]))
if (sum(r_inrange) == 0) {
  r_min <- -10
  r_max <- 10
}
r_vals <- seq(
  r_min,
  r_max,
  0.2
)
r_list <- lapply(numeric(2), function(x) r_vals)
r_mat <- as.matrix(expand.grid(r_list))

grid_mat <- expand_grid(time_vals, r_mat) %>%
  as.matrix()

sim_pred <- matrix(nrow = nrow(grid_mat), ncol = ncol(grid_mat))
for (i in 1:nrow(sim_pred)) {
  sim_pred[i, ] <- sim_result$embedding_map(unlist(as.vector(grid_mat[i, ])))
  progress(i, nrow(sim_pred))
}

sim_pred_full <- cbind(grid_mat, sim_pred)
sim_pred_full_df <- data.frame(sim_pred_full)
names(sim_pred_full_df) <- c("time", "r1", "r2", "x1", "x2", "x3")

range_x1 <- max(sim_pred_full_df$x1) - min(sim_pred_full_df$x1)
range_x2 <- max(sim_pred_full_df$x2) - min(sim_pred_full_df$x2)
range_x3 <- max(sim_pred_full_df$x3) - min(sim_pred_full_df$x3)

plot_ly(
  sim_pred_full_df,
  x = ~x1,
  y = ~x2,
  z = ~x3,
  frame = ~time,
  opacity = 0.5,
  type = "scatter3d",
  mode = "markers"
) %>%
  layout(
    scene = list(
      xaxis = list(
        range = c(
          min(sim_pred_full_df$x1) - (0.1 * range_x1),
          max(sim_pred_full_df$x1) + (0.1 * range_x1)
        )
      ),
      yaxis = list(
        range = c(
          min(sim_pred_full_df$x2) - (0.1 * range_x2),
          max(sim_pred_full_df$x2) + (0.1 * range_x2)
        )
      ),
      zaxis = list(
        range = c(
          min(sim_pred_full_df$x3) - (0.1 * range_x3),
          max(sim_pred_full_df$x3) + (0.1 * range_x3)
        )
      )
    )
  )
