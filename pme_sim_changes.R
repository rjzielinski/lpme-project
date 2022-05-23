library(gridExtra)
library(plotly)
library(tidyverse)

source("Principal_Manifold_Estimation.R")

sim_D2d1_case1 <- function(time_val, vertical_multiplier, horizontal_multiplier, noise, time_noise) { I <- 1000
  t <- rnorm(I, mean = 0, sd = 1)
  horizontal_noise <- rnorm(1, mean = 0, sd = time_noise)
  vertical_noise <- rnorm(1, mean = 0, sd = time_noise)
  sd.noise <- noise
  e1 <- rnorm(I, mean = 0, sd = sd.noise)
  e2 <- rnorm(I, mean = 0, sd = sd.noise)
  X <- matrix(NA, nrow = I, ncol = 2)
  manifold <- function(tau, time_val, vertical_multiplier, horizontal_multiplier, vertical_noise, horizontal_noise) {
    return(
      c(
        # tau + (horizontal_multiplier * sin(time_val)),
        tau + horizontal_noise,
        # sin(tau + pi / 2) + (vertical_multiplier * sin(time_val))
        sin(tau + pi / 2) + vertical_noise
      )
    )
  }
  
  X <- map(
    t, 
    ~ manifold(
      .x, 
      time_val, 
      vertical_multiplier, 
      horizontal_multiplier, 
      vertical_noise, 
      horizontal_noise
    )
  ) %>% 
    unlist() %>% 
    matrix(ncol = 2, byrow = TRUE)
  data.points <- X + cbind(e1, e2)
  
  result <- PME(x.obs = data.points, d = 1)
  
  f <- result$embedding.map
  t.test <- seq(from = -100, to = 100, by = 0.05)
  t.length <- length(t.test)
  x.test <- map(t.test, ~ f(.x)) %>% 
    unlist() %>% 
    matrix(ncol = 2, byrow = TRUE)
  
  p <- ggplot() +
    geom_point(aes(x = data.points[, 1], y = data.points[, 2]), color = "grey") +
    geom_line(aes(x = x.test[, 1], y = x.test[, 2]), color = "red") +
    xlim(min(data.points[, 1]) * 1.5, max(data.points[, 1]) * 1.5) +
    ylim(min(data.points[, 2]) * 1.5, max(data.points[, 2]) * 1.5)
  
  return(list(result, data.points, x.test, p))
}

time_vals <- seq(from = 0, to = 5, by = 0.5)

sim_D2d1_case1_results <- lapply(
  time_vals, 
  sim_D2d1_case1, 
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  noise = 0.15
)

saveRDS(sim_D2d1_case1_results, file = "sim_D2d1_case1_results.RDS")

sim_D2d1_case1_results_noise <- lapply(
  time_vals,
  sim_D2d1_case1,
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  noise = 0.75
)

time_val_vec <- rep(time_vals, each = dim(sim_D2d1_case1_results[[1]][[3]])[1])

sim_D2d1_case1_data_df <- map(
  1:length(time_vals),
  ~ data.frame(cbind(time_vals[.x], sim_D2d1_case1_results[[.x]][[2]]))
) %>% 
  reduce(bind_rows)
names(sim_D2d1_case1_data_df) <- c(
  "time",
  "x",
  "y"
)

sim_D2d1_case1_data_noise_df <- map(
  1:length(time_vals),
  ~ data.frame(cbind(time_vals[.x], sim_D2d1_case1_results_noise[[.x]][[2]]))
) %>% 
  reduce(bind_rows)
names(sim_D2d1_case1_data_noise_df) <- c(
  "time",
  "x",
  "y"
)

sim_D2d1_case1_min_x <- min(sim_D2d1_case1_data_df$x)
sim_D2d1_case1_max_x <- max(sim_D2d1_case1_data_df$x)
sim_D2d1_case1_range_x <- sim_D2d1_case1_max_x - sim_D2d1_case1_min_x
sim_D2d1_case1_min_y <- min(sim_D2d1_case1_data_df$y)
sim_D2d1_case1_max_y <- max(sim_D2d1_case1_data_df$y)
sim_D2d1_case1_range_y <- sim_D2d1_case1_max_y - sim_D2d1_case1_min_y

sim_D2d1_case1_noise_min_x <- min(sim_D2d1_case1_data_noise_df$x)
sim_D2d1_case1_noise_max_x <- max(sim_D2d1_case1_data_noise_df$x)
sim_D2d1_case1_noise_range_x <- sim_D2d1_case1_noise_max_x - sim_D2d1_case1_noise_min_x
sim_D2d1_case1_noise_min_y <- min(sim_D2d1_case1_data_noise_df$y)
sim_D2d1_case1_noise_max_y <- max(sim_D2d1_case1_data_noise_df$y)
sim_D2d1_case1_noise_range_y <- sim_D2d1_case1_noise_max_y - sim_D2d1_case1_noise_min_y


sim_D2d1_case1_manifold_df <- map(
  1:length(time_vals),
  ~ data.frame(cbind(time_vals[.x], sim_D2d1_case1_results[[.x]][[3]]))
) %>% 
  reduce(bind_rows)
names(sim_D2d1_case1_manifold_df) <- c(
  "time",
  "x",
  "y"
)

sim_D2d1_case1_noise_manifold_df <- map(
  1:length(time_vals),
  ~ data.frame(cbind(time_vals[.x], sim_D2d1_case1_results_noise[[.x]][[3]]))
) %>% 
  reduce(bind_rows)
names(sim_D2d1_case1_noise_manifold_df) <- c(
  "time",
  "x",
  "y"
)

sim_D2d1_case1_manifold_df_red <- sim_D2d1_case1_manifold_df %>% 
  filter(
    x > sim_D2d1_case1_min_x - (0.25 * sim_D2d1_case1_range_x),
    x < sim_D2d1_case1_max_x + (0.25 * sim_D2d1_case1_range_x),
    y > sim_D2d1_case1_min_y - (0.25 * sim_D2d1_case1_range_y),
    y < sim_D2d1_case1_max_y + (0.25 * sim_D2d1_case1_range_y)
  )

sim_D2d1_case1_noise_manifold_df_red <- sim_D2d1_case1_noise_manifold_df %>% 
  filter(
    x > sim_D2d1_case1_noise_min_x - (0.25 * sim_D2d1_case1_noise_range_x),
    x < sim_D2d1_case1_noise_max_x + (0.25 * sim_D2d1_case1_noise_range_x),
    y > sim_D2d1_case1_noise_min_y - (0.25 * sim_D2d1_case1_noise_range_y),
    y < sim_D2d1_case1_noise_max_y + (0.25 * sim_D2d1_case1_noise_range_y)
  )

plot_D2d1_case1 <- plot_ly(
  sim_D2d1_case1_manifold_df_red,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "lines",
  color = ~time
)

save_image(plot_D2d1_case1, "sim_D2d1_case1_line_plot.png")

plot_D2d1_case1_noise <- plot_ly(
  sim_D2d1_case1_noise_manifold_df_red,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "lines",
  color = ~time
)

sim_D2d1_case2 <- function(time_val, vertical_multiplier, horizontal_multiplier, noise) {
  manifold <- function(tau, time_val, vertical_multiplier, horizontal_multiplier) {
    return(
      c(
        cos(tau) + (horizontal_multiplier * sin(time_val)),
        sin(tau) + (vertical_multiplier * sin(time_val))
      )
    )
  }
  I <- 1000
  t <- runif(I, min = 0, max = 1.5 * pi)
  X <- manifold(t, time_val, vertical_multiplier, horizontal_multiplier)
  sd.noise <- noise
  e1 <- rnorm(I, mean = 0, sd = sd.noise)
  e2 <- rnorm(I, mean = 0, sd = sd.noise)
  
  data.points <- X + cbind(e1, e2)
  
  result <- PME(x.obs = data.points, d = 1)
  f <- result$embedding.map
  
  t.test <- seq(from = -100, to = 100, by = 0.05)
  t.length <- length(t.test)
  x.test <- map(t.test, ~ f(.x)) %>% 
    unlist() %>% 
    matrix(ncol = 2, byrow = TRUE)
  
  p <- ggplot() +
    geom_point(aes(x = data.points[, 1], y = data.points[, 2]), color = "grey") +
    geom_line(aes(x = x.test[, 1], y = x.test[, 2]), color = "red") +
    xlim(min(data.points[, 1]) * 1.5, max(data.points[, 1]) * 1.5) +
    ylim(min(data.points[, 2]) * 1.5, max(data.points[, 2]) * 1.5)
  
  return(list(result, data.points, x.test, p))
}

time_vals <- seq(from = 0, to = 5, by = 0.5)

sim_D2d1_case2_results <- lapply(
  time_vals, 
  sim_D2d1_case2, 
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  noise = 0.1
)

saveRDS(sim_D2d1_case2_results, "sim_D2d1_case2_results.RDS")

sim_D2d1_case2_results_noise <- lapply(
  time_vals,
  sim_D2d1_case2,
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  noise = 0.5
)

time_val_vec <- rep(time_vals, each = dim(sim_D2d1_case2_results[[1]][[3]])[1])

sim_D2d1_case2_data_df <- map(
  1:length(time_vals),
  ~ data.frame(cbind(time_vals[.x], sim_D2d1_case2_results[[.x]][[2]]))
) %>% 
  reduce(bind_rows)
names(sim_D2d1_case2_data_df) <- c(
  "time",
  "x",
  "y"
)

sim_D2d1_case2_data_noise_df <- map(
  1:length(time_vals),
  ~ data.frame(cbind(time_vals[.x], sim_D2d1_case2_results_noise[[.x]][[2]]))
) %>% 
  reduce(bind_rows)
names(sim_D2d1_case2_data_noise_df) <- c(
  "time",
  "x",
  "y"
)

sim_D2d1_case2_min_x <- min(sim_D2d1_case2_data_df$x)
sim_D2d1_case2_max_x <- max(sim_D2d1_case2_data_df$x)
sim_D2d1_case2_range_x <- sim_D2d1_case2_max_x - sim_D2d1_case2_min_x
sim_D2d1_case2_min_y <- min(sim_D2d1_case2_data_df$y)
sim_D2d1_case2_max_y <- max(sim_D2d1_case2_data_df$y)
sim_D2d1_case2_range_y <- sim_D2d1_case2_max_y - sim_D2d1_case2_min_y

sim_D2d1_case2_noise_min_x <- min(sim_D2d1_case2_data_noise_df$x)
sim_D2d1_case2_noise_max_x <- max(sim_D2d1_case2_data_noise_df$x)
sim_D2d1_case2_noise_range_x <- sim_D2d1_case2_noise_max_x - sim_D2d1_case2_noise_min_x
sim_D2d1_case2_noise_min_y <- min(sim_D2d1_case2_data_noise_df$y)
sim_D2d1_case2_noise_max_y <- max(sim_D2d1_case2_data_noise_df$y)
sim_D2d1_case2_noise_range_y <- sim_D2d1_case2_noise_max_y - sim_D2d1_case2_noise_min_y


sim_D2d1_case2_manifold_df <- map(
  1:length(time_vals),
  ~ data.frame(cbind(time_vals[.x], sim_D2d1_case2_results[[.x]][[3]]))
) %>% 
  reduce(bind_rows)
names(sim_D2d1_case2_manifold_df) <- c(
  "time",
  "x",
  "y"
)

sim_D2d1_case2_noise_manifold_df <- map(
  1:length(time_vals),
  ~ data.frame(cbind(time_vals[.x], sim_D2d1_case2_results_noise[[.x]][[3]]))
) %>% 
  reduce(bind_rows)
names(sim_D2d1_case2_noise_manifold_df) <- c(
  "time",
  "x",
  "y"
)

sim_D2d1_case2_manifold_df_red <- sim_D2d1_case2_manifold_df %>% 
  filter(
    x > sim_D2d1_case2_min_x - (0.25 * sim_D2d1_case2_range_x),
    x < sim_D2d1_case2_max_x + (0.25 * sim_D2d1_case2_range_x),
    y > sim_D2d1_case2_min_y - (0.25 * sim_D2d1_case2_range_y),
    y < sim_D2d1_case2_max_y + (0.25 * sim_D2d1_case2_range_y)
  )

sim_D2d1_case2_noise_manifold_df_red <- sim_D2d1_case2_noise_manifold_df %>% 
  filter(
    x > sim_D2d1_case2_noise_min_x - (0.25 * sim_D2d1_case2_noise_range_x),
    x < sim_D2d1_case2_noise_max_x + (0.25 * sim_D2d1_case2_noise_range_x),
    y > sim_D2d1_case2_noise_min_y - (0.25 * sim_D2d1_case2_noise_range_y),
    y < sim_D2d1_case2_noise_max_y + (0.25 * sim_D2d1_case2_noise_range_y)
  )

plot_D2d1_case2 <- plot_ly(
  sim_D2d1_case2_manifold_df_red,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "lines",
  color = ~time
)

save_image(plot_D2d1_case2, "sim_D2d1_case2_line_plot.png")

plot_D2d1_case2_noise <- plot_ly(
  sim_D2d1_case2_noise_manifold_df_red,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "lines",
  color = ~time
)

### Case I

sim_D3d1_case1 <- function(time_val, vertical_multiplier, horizontal_multiplier, noise) {
  manifold <- function(tau, time_val, vertical_multiplier, horizontal_multiplier) {
    return(
      c(
        tau + time_val, 
        tau^2 + (horizontal_multiplier * sin(time_val)),
        tau^3 + (vertical_multiplier * sin(time_val))
      )
    )
  }
  
  t <- runif(I, min = -1, max = 1)
  
  X <- manifold(t, time_val, vertical_multiplier, horizontal_multiplier)
  e1 <- rnorm(I, mean = 0, sd = noise)
  e2 <- rnorm(I, mean = 0, sd = noise)
  e3 <- rnorm(I, mean = 0, sd = noise)
  data.points <- X + cbind(e1, e2, e3)
  
  result <- PME(x.obs = data.points, d = 1)
  f <- result$embedding.map
  
  t.test <- seq(from = -100, to = 100, by = 0.05)
  t.length <- length(t.test)
  x.test <- map(t.test, ~ f(.x)) %>% 
    unlist() %>% 
    matrix(ncol = 3, byrow = TRUE)

f=result$embedding.map
t.test=seq(from=-10,to=10,by=0.005)
t.length=length(t.test)
x.test=matrix(NA,ncol=3,nrow = t.length)
for(i in 1:t.length){ x.test[i,]=f(t.test[i]) }
index=(x.test[,1]>=min(data.points[,1]))&x.test[,1]<=max(data.points[,1])&(x.test[,2]>=min(data.points[,2]))&x.test[,2]<=max(data.points[,2])&(x.test[,3]>=min(data.points[,3]))&x.test[,3]<=max(data.points[,3])
scatter3D(data.points[,1], data.points[,2], data.points[,3], 
          pch = 20, box=TRUE, cex = 0.2, colkey = FALSE, 
          border="black", shade=0.8, 
          bty = "g", ticktype = "detailed",
          main="Principal Manifold Estimation")
scatter3D(x.test[index,1], x.test[index,2],x.test[index,3], 
          pch = 20, box=TRUE, cex = 0.2, colkey = FALSE, col = "red", 
          border="black", shade=0.8, main=" ",add = TRUE)

  result <- PME(x.obs = data.points, d = 1)
  
  
  
  p <- ggplot() +
    geom_point(aes(x = data.points[, 1], y = data.points[, 2]), color = "grey") +
    geom_line(aes(x = x.test[, 1], y = x.test[, 2]), color = "red") +
    xlim(min(data.points[, 1]) * 1.5, max(data.points[, 1]) * 1.5) +
    ylim(min(data.points[, 2]) * 1.5, max(data.points[, 2]) * 1.5)
  
  return(list(result, data.points, x.test, p))
}

time_vals <- seq(from = 0, to = 5, by = 0.5)

sim_D2d1_case2_results <- lapply(
  time_vals, 
  sim_D2d1_case2, 
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  noise = 0.1
)

sim_D2d1_case2_results_noise <- lapply(
  time_vals,
  sim_D2d1_case2,
  vertical_multiplier = 1,
  horizontal_multiplier = 1,
  noise = 0.5
)

time_val_vec <- rep(time_vals, each = dim(sim_D2d1_case2_results[[1]][[3]])[1])

sim_D2d1_case2_data_df <- map(
  1:length(time_vals),
  ~ data.frame(cbind(time_vals[.x], sim_D2d1_case2_results[[.x]][[2]]))
) %>% 
  reduce(bind_rows)
names(sim_D2d1_case2_data_df) <- c(
  "time",
  "x",
  "y"
)

sim_D2d1_case2_data_noise_df <- map(
  1:length(time_vals),
  ~ data.frame(cbind(time_vals[.x], sim_D2d1_case2_results_noise[[.x]][[2]]))
) %>% 
  reduce(bind_rows)
names(sim_D2d1_case2_data_noise_df) <- c(
  "time",
  "x",
  "y"
)

sim_D2d1_case2_min_x <- min(sim_D2d1_case2_data_df$x)
sim_D2d1_case2_max_x <- max(sim_D2d1_case2_data_df$x)
sim_D2d1_case2_range_x <- sim_D2d1_case2_max_x - sim_D2d1_case2_min_x
sim_D2d1_case2_min_y <- min(sim_D2d1_case2_data_df$y)
sim_D2d1_case2_max_y <- max(sim_D2d1_case2_data_df$y)
sim_D2d1_case2_range_y <- sim_D2d1_case2_max_y - sim_D2d1_case2_min_y

sim_D2d1_case2_noise_min_x <- min(sim_D2d1_case2_data_noise_df$x)
sim_D2d1_case2_noise_max_x <- max(sim_D2d1_case2_data_noise_df$x)
sim_D2d1_case2_noise_range_x <- sim_D2d1_case2_noise_max_x - sim_D2d1_case2_noise_min_x
sim_D2d1_case2_noise_min_y <- min(sim_D2d1_case2_data_noise_df$y)
sim_D2d1_case2_noise_max_y <- max(sim_D2d1_case2_data_noise_df$y)
sim_D2d1_case2_noise_range_y <- sim_D2d1_case2_noise_max_y - sim_D2d1_case2_noise_min_y


sim_D2d1_case2_manifold_df <- map(
  1:length(time_vals),
  ~ data.frame(cbind(time_vals[.x], sim_D2d1_case2_results[[.x]][[3]]))
) %>% 
  reduce(bind_rows)
names(sim_D2d1_case2_manifold_df) <- c(
  "time",
  "x",
  "y"
)

sim_D2d1_case2_noise_manifold_df <- map(
  1:length(time_vals),
  ~ data.frame(cbind(time_vals[.x], sim_D2d1_case2_results_noise[[.x]][[3]]))
) %>% 
  reduce(bind_rows)
names(sim_D2d1_case2_noise_manifold_df) <- c(
  "time",
  "x",
  "y"
)

sim_D2d1_case2_manifold_df_red <- sim_D2d1_case2_manifold_df %>% 
  filter(
    x > sim_D2d1_case2_min_x - (0.25 * sim_D2d1_case2_range_x),
    x < sim_D2d1_case2_max_x + (0.25 * sim_D2d1_case2_range_x),
    y > sim_D2d1_case2_min_y - (0.25 * sim_D2d1_case2_range_y),
    y < sim_D2d1_case2_max_y + (0.25 * sim_D2d1_case2_range_y)
  )

sim_D2d1_case2_noise_manifold_df_red <- sim_D2d1_case2_noise_manifold_df %>% 
  filter(
    x > sim_D2d1_case2_noise_min_x - (0.25 * sim_D2d1_case2_noise_range_x),
    x < sim_D2d1_case2_noise_max_x + (0.25 * sim_D2d1_case2_noise_range_x),
    y > sim_D2d1_case2_noise_min_y - (0.25 * sim_D2d1_case2_noise_range_y),
    y < sim_D2d1_case2_noise_max_y + (0.25 * sim_D2d1_case2_noise_range_y)
  )

plot_D2d1_case2 <- plot_ly(
  sim_D2d1_case2_manifold_df_red,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "lines",
  color = ~time
)

plot_D2d1_case2_noise <- plot_ly(
  sim_D2d1_case2_noise_manifold_df_red,
  x = ~x,
  y = ~y,
  z = ~time,
  type = "scatter3d",
  mode = "lines",
  color = ~time
)



### Case II

I=1000
manifold=function(t){ return(c(t,cos(t),sin(t))) }
t=seq(from=0,to=3*pi,length.out = I)
X=matrix(0,nrow = length(t),ncol = 3)
for(i in 1:length(t)){ X[i,]=manifold(t[i]) }
noise=0.05
e1=rnorm(I,mean=0,sd=noise)
e2=rnorm(I,mean=0,sd=noise)
e3=rnorm(I,mean=0,sd=noise)
data.points=X+cbind(e1,e2,e3)

ptm <- proc.time()
result=PME(x.obs=data.points, d=1)
proc.time() - ptm

f=result$embedding.map
t.test=seq(from=-20,to=20,by=0.01)
t.length=length(t.test)
x.test=matrix(NA,ncol=3,nrow = t.length)
for(i in 1:t.length){ x.test[i,]=f(t.test[i]) }
index=(x.test[,1]>=min(data.points[,1]))&x.test[,1]<=max(data.points[,1])&(x.test[,2]>=min(data.points[,2]))&x.test[,2]<=max(data.points[,2])&(x.test[,3]>=min(data.points[,3]))&x.test[,3]<=max(data.points[,3])
scatter3D(data.points[,1], data.points[,2], data.points[,3], 
          pch = 20, box=TRUE, cex = 0.2, colkey = FALSE, 
          border="black", shade=0.8, 
          bty = "g", ticktype = "detailed",
          main="Principal Manifold Estimation")
scatter3D(x.test[index,1], x.test[index,2],x.test[index,3], 
          pch = 20, box=TRUE, cex = 0.2, colkey = FALSE, col = "red", 
          border="black", shade=0.8, main=" ",add = TRUE)
