library(tidyverse)
library(plotly)
library(pracma)
library(profvis)

source("code/lpme.R")
source("code/pme.R")
source("code/functions/lpme_init.R")

sim_data <- function(time_val, case, noise, shape_noise) {
  manifolds <- list(
    function(tau, amp_noise, period_noise) {
      return(c(tau, amp_noise[1] * sin(period_noise[1] * tau + pi / 2)))
    },
    function(tau, amp_noise, period_noise) {
      return(c(tau, amp_noise[1] * sin(period_noise[1] * tau)))
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          amp_noise[1] * cos(period_noise[1] * tau),
          amp_noise[2] * sin(period_noise[2] * tau)
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          amp_noise[1] * cos(period_noise[1] * tau),
          amp_noise[2] * sin(period_noise[2] * tau)
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          tau,
          (tau * amp_noise[1] + period_noise[1]) ^ 2,
          (tau * amp_noise[2] + period_noise[2]) ^ 3
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          tau,
          amp_noise[1] * cos(period_noise[1] * tau),
          amp_noise[2] * sin(period_noise[2] * tau)
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          period_noise[1] * tau[1],
          period_noise[2] * tau[2],
          amp_noise[1] * (amp_noise[2] * norm_euclidean(period_noise * tau) ^ 2)
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          (period_noise[1] * amp_noise[1] * tau[1]) * cos(amp_noise[1] * tau[1]),
          (period_noise[2] * amp_noise[2] * tau[1]) * sin(amp_noise[2] * tau[1]),
          tau[2]
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          amp_noise[1] * sin(period_noise[1] * tau[1]) * cos(period_noise[2] * tau[2]),
          amp_noise[1] * sin(period_noise[1] * tau[1]) * sin(period_noise[2] * tau[2]),
          amp_noise[1] * cos(period_noise[1] * tau[1])
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      r <- 1 + amp_noise[1] * (tau[1] + 1) * sqrt(tau[2] + 1)
      return(
        c(
          r * sin(period_noise[1] * tau[1]) * cos(period_noise[2] * tau[2]),
          r * sin(period_noise[1] * tau[1]) * sin(period_noise[2] * tau[2]),
          r * cos(period_noise[1] * tau[1])
        )
      )
    }
  )

  D <- case_when(
    case == 1 ~ 2,
    case == 2 ~ 2,
    case == 3 ~ 2,
    case == 4 ~ 2,
    case == 5 ~ 3,
    case == 6 ~ 3,
    case == 7 ~ 3,
    case == 8 ~ 3,
    case == 9 ~ 3,
    case == 10 ~ 3
  )
  d <- case_when(
    case == 1 ~ 1,
    case == 2 ~ 1,
    case == 3 ~ 1,
    case == 4 ~ 1,
    case == 5 ~ 1,
    case == 6 ~ 1,
    case == 7 ~ 2,
    case == 8 ~ 2,
    case == 9 ~ 2,
    case == 10 ~ 2
  )

  manifold <- manifolds[[case]]

  I <- 1000
  t <- matrix(NA, nrow = I, ncol = d)
  X <- matrix(NA, nrow = I, ncol = D)
  noise_vals <- rnorm(I * D, mean = 0, sd = noise) %>%
    matrix(nrow = I, ncol = D)

  if (case == 1) {
    t[, 1] <- rnorm(I, mean = 0, sd = 1)
  } else if (case == 2) {
    t[, 1] <- runif(I, min = -3 * pi, max = 3 * pi)
  } else if (case == 3) {
    t[, 1] <- runif(I, min = 0, max = 1.5 * pi)
  } else if (case == 4) {
    t[, 1] <- runif(I, time_val * (pi / 4), (time_val + 6) * (pi / 4))
  } else if (case == 5) {
    t[, 1] <- runif(I, min = -1, max = 1)
  } else if (case == 6) {
    t[, 1] <- runif(I, min = 0, max = 3 * pi)
  } else if (case == 7) {
    t[, 1] <- runif(I, min = -1, max = 1)
    t[, 2] <- runif(I, min = -1, max = 1)
  } else if (case == 8) {
    t[, 1] <- runif(I, min = 0, max = 10)
    t[, 2] <- runif(I, min = -1, max = 1)
  } else if (case == 9) {
    t[, 1] <- runif(I, min = 0, max = pi)
    t[, 2] <- runif(I, min = 0, max = 2 * pi)
  } else if (case == 10) {
    t[, 1] <- runif(I, min = 0, max = pi)
    t[, 2] <- runif(I, min = 0, max = 2 * pi)
  }

  if (case %in% 1:10) {
    amp_mean <- 1
    per_mean <- 1
  }

  amp_noise <- rnorm(d, mean = amp_mean, sd = shape_noise)
  period_noise <- rnorm(d, mean = per_mean, sd = shape_noise)

  X <- map(
    1:nrow(t),
    ~ manifold(
      t[.x, ],
      amp_noise,
      period_noise
    )
  ) %>%
    unlist() %>%
    matrix(ncol = D, byrow = TRUE)

  data.points <- X + noise_vals
  data.points <- cbind(time_val, data.points)
  return(data.points)
}

########## INITIALIZATION SIMULATIONS ##########

time_noise_vals <- c(0.1, 0.25, 0.5, 1)
t_vals <- c(5, 10, 25)
functions <- list(
  sim_D2d1_case1,
  sim_D2d1_case2,
  sim_D2d1_case3,
  sim_D2d1_case4,
  sim_D3d1_case1,
  sim_D3d1_case2,
  sim_D3d2_case1,
  sim_D3d2_case2
)

dim_names <- c(
  "D2d1",
  "D2d1",
  "D2d1",
  "D2d1",
  "D3d1",
  "D3d1",
  "D3d2",
  "D3d2"
)

dim_vals <- c(1, 1, 1, 1, 1, 1, 2, 2)

case_names <- c(
  "case1",
  "case2",
  "case3",
  "case4",
  "case1",
  "case2",
  "case1",
  "case2"
)

init_types <- c("first", "full", "separate")

set.seed(500)
for (noise_val in time_noise_vals) {
  for (t_val in t_vals) {
    for (fun_idx in 1:length(functions)) {
      time_vals <- 0:t_val
      df_list <- lapply(
        time_vals,
        functions[[fun_idx]],
        vertical_multiplier = 0.1,
        horizontal_multiplier = 0.1,
        noise = 0.15,
        time_noise = noise_val
      )
      data_points <- reduce(df_list, rbind)
      for (init_type in init_types) {
        sim_init_result <- lpme_init(
          data_points,
          dim_vals[fun_idx],
          init = init_type
        )

        init_result_dir <- paste0(
          "results/init_sim/",
          dim_names[fun_idx],
          "/noise_",
          gsub(pattern = "\\.", replacement = "", x = as.character(noise_val)),
          "/t",
          as.character(length(time_vals)),
          "/",
          case_names[fun_idx],
          "_",
          init_type
        )

        if (!dir.exists(init_result_dir)) {
          dir.create(init_result_dir, recursive = TRUE)
        }

        saveRDS(
          sim_init_result,
          paste0(init_result_dir, ".RDS")
        )

        sim_lpme_result <- lpme(
          data_points,
          dim_vals[fun_idx],
          init = init_type
        )

        lpme_result_dir <- paste0(
          "results/lpme/",
          dim_names[fun_idx],
          "/noise_",
          gsub(pattern = "\\.", replacement = "", x = as.character(noise_val)),
          "/t",
          as.character(length(time_vals)),
          "/",
          case_names[fun_idx],
          "_",
          init_type
        )

        if (!dir.exists(lpme_result_dir)) {
          dir.create(lpme_result_dir, recursive = TRUE)
        }

        saveRDS(
          sim_lpme_result,
          paste0(lpme_result_dir, ".RDS")
        )
      }
    }
  }
}

### Simulation Case 1

time_vals <- seq(0, 10, 2)

set.seed(100)
df_list <- lapply(
  time_vals,
  sim_D2d1_case1,
  vertical_multiplier = 0,
  horizontal_multiplier = 0,
  noise = 0.15,
  time_noise = 0.1,
  shape_noise = 0.3
)

data_points <- reduce(df_list, rbind)

sim_result <- lpme(data_points, 1)

time_vals <- seq(0, 10, 0.1)
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
  mode = "markers",
  marker = list(
    size = 1
  )
) %>%
  add_trace(
    x = data_points[, 2],
    y = data_points[, 3],
    z = data_points[, 1]
  )

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

