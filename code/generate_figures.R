library(plot3D)
library(plotly)
library(pracma)
library(RColorBrewer)
library(scatterplot3d)
library(tidyverse)

set.seed(1202)
theta <- runif(1000, min = 0, max = 2 * pi)
x <- cos(theta)
y <- sin(theta)

unit_df <- data.frame(theta, x, y) %>%
  arrange(theta)
angle <- cart2pol(as.matrix(unit_df[, 2:3], ncol = 2))


png("paper/figures/unit_circle_augmentation.png", res = 500, height = 2000, width = 2000)
par(mfrow = c(1, 2))
plot(
  x = unit_df$x,
  y = unit_df$y,
  type = "l",
  asp = 1,
  main = "Unit Circle, D = 2",
  xlab = "x",
  ylab = "y"
)

plt2 <- lines3D(
  x = unit_df$x,
  y = unit_df$y,
  z = angle[, 1],
  main = "Unit Circle, D = 3",
  xlab = "x",
  ylab = "y",
  zlab = ""
)
text3D(-2, -0.8, -2.25, labels = expression(theta), add = TRUE)
dev.off()

png("paper/figures/unit_circle_D2.png", res = 500, height = 2000, width = 2000)
plot(
  x = unit_df$x,
  y = unit_df$y,
  type = "l",
  asp = 1,
  main = "Unit Circle, D = 2",
  xlab = "x",
  ylab = "y"
)
dev.off()

png("paper/figures/unit_circle_D3.png", res = 500, height = 2000, width = 2000)
lines3D(
  x = unit_df$x,
  y = unit_df$y,
  z = angle[, 1],
  main = "Unit Circle, D = 3",
  xlab = "x",
  ylab = "y",
  zlab = ""
)
text3D(-2.25, -0.85, -2.25, labels = expression(theta), add = TRUE)
dev.off()

sim_example_case1 <- readRDS("simulations/case1/duration_01_interval_010_ampnoise_025_pernoise_050_n_1000_constant000_run_01.rds")
sim_example_case5 <- readRDS("simulations/case6/duration_01_interval_010_ampnoise_000_pernoise_005_n_1000_constant000_run_01.rds")
sim_example_case7 <- readRDS("simulations/case7/duration_01_interval_010_ampnoise_010_pernoise_010_n_1000_constant000_run_01.rds")

df_labels <- c("Data", "True", "LPME", "PME", "Principal Curve")

data_case1 <- list(
  sim_example_case1$df,
  sim_example_case1$true_vals,
  sim_example_case1$lpme_vals,
  sim_example_case1$pme_vals,
  sim_example_case1$principal_curve_vals
)

for (df_idx in 1:length(data_case1)) {
  df <- as.data.frame(data_case1[[df_idx]])
  names(df) <- c("time", "x", "y")
  df$type <- df_labels[df_idx]
  data_case1[[df_idx]] <- df
}
df_case1 <- reduce(data_case1, rbind)
df_case1 <- df_case1 %>%
  filter(mod(time * 10, 2) == 0)

time_vals <- sim_example_case1$times
time_vals <- time_vals[mod(time_vals * 10, 2) == 0]

png("paper/figures/sim_case1.png", res = 500, height = 2000, width = 3000)
par(oma = c(4, 1, 1, 1), mfrow = c(2, 3), mar = c(2, 2, 1, 1))
for (time_idx in 1:length(time_vals)) {
  temp_data <- df_case1 %>%
    filter(
      time == time_vals[time_idx],
      type == "Data"
    ) %>%
    arrange(x)
  # temp_data <- temp_data[temp_data[, 1] == time_vals[time_idx], ]
  # temp_lpme <- sim_example_case1$lpme_vals[temp_data[, 1] == time_vals[time_idx], ]
  temp_lpme <- df_case1 %>%
    filter(
      time == time_vals[time_idx],
      type == "LPME"
    ) %>%
    arrange(x)
  temp_pme <- df_case1 %>%
    filter(
      time == time_vals[time_idx],
      type == "PME"
    ) %>%
    arrange(x)
  temp_princurve <- df_case1 %>%
    filter(
      time == time_vals[time_idx],
      type == "Principal Curve"
    ) %>%
    arrange(x)
  temp_true <- df_case1 %>%
    filter(
      time == time_vals[time_idx],
      type == "True"
    ) %>%
    arrange(x)
  # temp_pme <- sim_example_case1$pme_vals[temp_data[, 1] == time_vals[time_idx], ]
  # temp_princurve <- sim_example_case1$principal_curve_vals[temp_data[, 1] == time_vals[time_idx], ]
  # temp_true <- sim_example_case1$true_vals[temp_data[, 1] == time_vals[time_idx], ]
  plot(
    x = temp_data$x,
    y = temp_data$y,
    xlab = "x",
    ylab = "y",
    col = alpha("black", 0.2),
    bg = alpha("black", 0.2),
    xlim = c(-4, 4),
    ylim = c(-2, 2),
    pch = 21,
    main = paste0("Time = ", time_vals[time_idx]),
    cex = 0.5
  )
  points(
    x = temp_true$x,
    y = temp_true$y,
    col = "red",
    bg = "red",
    pch = 21,
    # cex = 0.5
    type = "l",
    lwd = 1.25
  )
  points(
    x = temp_lpme$x,
    y = temp_lpme$y,
    col = "blue",
    bg = "blue",
    pch = 21,
    type = "l",
    lty = "dashed",
    lwd = 1.25
  )
  points(
    x = temp_pme$x,
    y = temp_pme$y,
    col = "green",
    bg = "green",
    pch = 21,
    type = "l",
    lty = "dotted",
    lwd = 1.25
  )
  points(
    x = temp_princurve$x,
    y = temp_princurve$y,
    col = "purple",
    bg = "purple",
    pch = 21,
    type = "l",
    lty = "dotdash",
    lwd = 1.25
  )
}
par(fig = c(0, 1, 0, 1), oma = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "l", bty = "n", xaxt = "n", yaxt = "n")
legend(
  "bottom",
  c("Data", "True Manifold", "LPME", "PME", "Principal Curve"),
  col = c("black", "red", "blue", "green", "purple"),
  # fill = c("black", "red", "blue", "green", "purple"),
  lty = c("solid", "solid", "dashed", "dotted", "dotdash"),
  lwd = 1.25,
  xpd = TRUE,
  horiz = TRUE,
  cex = 1,
  bty = "n"
)
dev.off()

data_case5 <- list(
  sim_example_case5$df,
  sim_example_case5$true_vals,
  sim_example_case5$lpme_vals,
  sim_example_case5$pme_vals,
  sim_example_case5$principal_curve_vals
)

for (df_idx in 1:length(data_case5)) {
  df <- as.data.frame(data_case5[[df_idx]])
  names(df) <- c("time", "x", "y", "z")
  df$type <- df_labels[df_idx]
  data_case5[[df_idx]] <- df
}
df_case5 <- reduce(data_case5, rbind)
df_case5 <- df_case5 %>%
  filter(mod(time * 10, 2) == 0)

time_vals <- sim_example_case5$times
time_vals <- time_vals[mod(time_vals * 10, 2) == 0]

png("paper/figures/sim_case5.png", res = 500, height = 2000, width = 3000)
par(oma = c(4, 1, 1, 1), mfrow = c(2, 3), mar = c(2, 2, 1, 1))
for (time_idx in 1:length(time_vals)) {
  temp_df <- df_case5 %>%
    filter(time == time_vals[time_idx]) %>%
    arrange(x)
  temp_data <- df_case5 %>%
    filter(
      time == time_vals[time_idx],
      type == "Data"
    ) %>%
    arrange(x)
  # temp_data <- temp_data[temp_data[, 1] == time_vals[time_idx], ]
  # temp_lpme <- sim_example_case1$lpme_vals[temp_data[, 1] == time_vals[time_idx], ]
  temp_lpme <- df_case5 %>%
    filter(
      time == time_vals[time_idx],
      type == "LPME"
    ) %>%
    arrange(x)
  temp_pme <- df_case5 %>%
    filter(
      time == time_vals[time_idx],
      type == "PME"
    ) %>%
    arrange(x)
  temp_princurve <- df_case5 %>%
    filter(
      time == time_vals[time_idx],
      type == "Principal Curve"
    ) %>%
    arrange(x)
  temp_true <- df_case5 %>%
    filter(
      time == time_vals[time_idx],
      type == "True"
    ) %>%
    arrange(x)
  # temp_pme <- sim_example_case1$pme_vals[temp_data[, 1] == time_vals[time_idx], ]
  # temp_princurve <- sim_example_case1$principal_curve_vals[temp_data[, 1] == time_vals[time_idx], ]
  # temp_true <- sim_example_case1$true_vals[temp_data[, 1] == time_vals[time_idx], ]

  plt <- scatter3D(
    x = temp_data$x,
    y = temp_data$y,
    z = temp_data$z,
    xlab = "x",
    ylab = "y",
    zlab = "z",
    xlim = c(-1, 10),
    ylim = c(-3, 3),
    zlim = c(-3, 3),
    pch = 21,
    col = alpha("black", 0.005),
    bg = alpha("black", 0.005),
    main = paste0("Time = ", time_vals[time_idx]),
    ticktype = "detailed"
  )
  scatter3D(
    x = temp_true$x,
    y = temp_true$y,
    z = temp_true$z,
    col = alpha("red", 0.5),
    bg = alpha("red", 0.5),
    pch = 20,
    type = "l",
    lwd = 1,
    add = TRUE
  )
  scatter3D(
    x = temp_lpme$x,
    y = temp_lpme$y,
    z = temp_lpme$z,
    col = "blue",
    bg = "blue",
    pch = 20,
    type = "l",
    lty = 2,
    lwd = 1,
    add = TRUE
  )
  scatter3D(
    x = temp_pme$x,
    y = temp_pme$y,
    z = temp_pme$z,
    col = "green",
    bg = "green",
    pch = 20,
    type = "l",
    lty = 3,
    lwd = 1,
    add = TRUE
  )
  scatter3D(
    x = temp_princurve$x,
    y = temp_princurve$y,
    z = temp_princurve$z,
    col = "purple",
    bg = "purple",
    pch = 20,
    type = "l",
    lty = 4,
    lwd = 1,
    add = TRUE
  )
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "l", bty = "n", xaxt = "n", yaxt = "n")
legend(
  "bottom",
  c("Data", "True Manifold", "LPME", "PME", "Principal Curve"),
  col = c("black", "red", "blue", "green", "purple"),
  # fill = c("black", "red", "blue", "green", "purple"),
  lwd = 1,
  xpd = TRUE,
  horiz = TRUE,
  cex = 1,
  bty = "n"
)
dev.off()

data_case7 <- list(
  sim_example_case7$df,
  sim_example_case7$true_vals,
  sim_example_case7$lpme_vals,
  sim_example_case7$pme_vals,
  sim_example_case7$principal_curve_vals
)

for (df_idx in 1:length(data_case7)) {
  df <- as.data.frame(data_case7[[df_idx]])
  if (dim(df)[2] == 4) {
    names(df) <- c("time", "x", "y", "z")
  } else if (dim(df)[2] == 5) {
    names(df) <- c("time", "x", "y", "z", "r")
  }
  df$type <- df_labels[df_idx]
  if ("r" %in% names(df)) {
    data_case7[[df_idx]] <- dplyr::select(df, -r)
  } else {
    data_case7[[df_idx]] <- df
  }
}
df_case7 <- reduce(data_case7, rbind)
df_case7 <- df_case7 %>%
  filter(mod(time * 10, 2) == 0)

time_vals <- unique(df_case7$time)

png("paper/figures/sim_case7.png", res = 500, height = 2000, width = 3000)
par(oma = c(4, 1, 1, 1), mfrow = c(2, 3), mar = c(2, 2, 1, 1))
for (time_idx in 1:length(time_vals)) {
  temp_data <- df_case7 %>%
    filter(
      time == time_vals[time_idx],
      type == "Data"
    ) %>%
    arrange(x)
  # temp_data <- temp_data[temp_data[, 1] == time_vals[time_idx], ]
  # temp_lpme <- sim_example_case1$lpme_vals[temp_data[, 1] == time_vals[time_idx], ]
  temp_lpme <- df_case7 %>%
    filter(
      time == time_vals[time_idx],
      type == "LPME"
    ) %>%
    arrange(x)
  temp_pme <- df_case7 %>%
    filter(
      time == time_vals[time_idx],
      type == "PME"
    ) %>%
    arrange(x)
  temp_princurve <- df_case7 %>%
    filter(
      time == time_vals[time_idx],
      type == "Principal Curve"
    ) %>%
    arrange(x)
  temp_true <- df_case7 %>%
    filter(
      time == time_vals[time_idx],
      type == "True"
    ) %>%
    arrange(x)
  # temp_pme <- sim_example_case1$pme_vals[temp_data[, 1] == time_vals[time_idx], ]
  # temp_princurve <- sim_example_case1$principal_curve_vals[temp_data[, 1] == time_vals[time_idx], ]
  # temp_true <- sim_example_case1$true_vals[temp_data[, 1] == time_vals[time_idx], ]

  plt <- scatter3D(
    x = temp_data$x,
    y = temp_data$y,
    z = temp_data$z,
    xlab = "x",
    ylab = "y",
    zlab = "z",
    xlim = c(-1, 1),
    ylim = c(-1, 1),
    zlim = c(-1, 3),
    pch = 21,
    col = alpha("black", 0.05),
    bg = alpha("black", 0.05),
    cex = 0.1,
    main = paste0("Time = ", time_vals[time_idx])
  )
  scatter3D(
    x = temp_true$x,
    y = temp_true$y,
    z = temp_true$z,
    col = alpha("red", 0.25),
    bg = alpha("red", 0.25),
    pch = 21,
    cex = 0.1,
    add = TRUE
  )
  scatter3D(
    x = temp_lpme$x,
    y = temp_lpme$y,
    z = temp_lpme$z,
    col = alpha("blue", 0.25),
    bg = alpha("blue", 0.25),
    pch = 21,
    cex = 0.1,
    add = TRUE
  )
  scatter3D(
    x = temp_pme$x,
    y = temp_pme$y,
    z = temp_pme$z,
    col = alpha("green", 0.25),
    bg = alpha("green", 0.25),
    pch = 21,
    cex = 0.1,
    add = TRUE
  )
  scatter3D(
    x = temp_princurve$x,
    y = temp_princurve$y,
    z = temp_princurve$z,
    col = alpha("purple", 0.25),
    bg = alpha("purple", 0.25),
    pch = 21,
    cex = 0.1,
    add = TRUE
  )
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "l", bty = "n", xaxt = "n", yaxt = "n")
legend(
  "bottom",
  c("Data", "True Manifold", "LPME", "PME", "Principal Curve"),
  col = c("black", "red", "blue", "green", "purple"),
  # fill = c("black", "red", "blue", "green", "purple"),
  lwd = 5,
  xpd = TRUE,
  horiz = TRUE,
  cex = 1,
  bty = "n"
)
dev.off()

### Following code follows from the sim_error_case1() function in lpme_simulations.R

set.seed(71927)
source("code/functions/sim_data.R")
max_time <- 3
interval <- 0.1
amp_noise <- 0.5
shape_noise <- 0.25
n <- 1000
time_change <- 0
time_trend <- "constant"

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

sim_df <- scale(
  sim_df,
  center = FALSE,
  scale = apply(sim_df, 2, function(x) max(abs(x)))
)
true_vals <- scale(
  true_vals,
  center = FALSE,
  scale = apply(sim_df, 2, function(x) max(abs(x)))
)
time_vals <- unique(sim_df[, 1])

set.seed(81752)
pme_result <- list()
pme_vals <- list()
for (t in 1:length(time_vals)) {
  print(t)
  temp_data <- sim_df[sim_df[, 1] == time_vals[t], -1]
  pme_result[[t]] <- pme(temp_data, d = 1, verbose = FALSE)
}

min_param <- min(
  map(
    1:length(time_vals),
    ~ min(pme_result[[.x]]$parameterization[[which.min(pme_result[[.x]]$MSD)]])
  ) %>%
    reduce(c)
)

max_param <- max(
  map(
    1:length(time_vals),
    ~ max(pme_result[[.x]]$parameterization[[which.min(pme_result[[.x]]$MSD)]])
  ) %>%
    reduce(c)
)

param_vals <- seq(min_param, max_param, length.out = 1000)

output_vals <- list()
for (i in 1:length(time_vals)) {
  output_vals[[i]] <- map(param_vals, ~ pme_result[[i]]$embedding_map(.x)) %>%
    reduce(rbind) %>%
    cbind(param_vals) %>%
    cbind(time_vals[i])
}
output_vals <- reduce(output_vals, rbind)

ggplot() +
  geom_line(
    aes(
      x = output_vals[, 1],
      y = output_vals[, 2],
      color = output_vals[, 3],
      group = output_vals[, 4]
    ),
    lwd = 2
  )

col_palette <- brewer.pal()
scatter3D(
  x = output_vals[, 1],
  z = output_vals[, 2],
  y = output_vals[, 3],
  colvar = output_vals[, 4],
  colkey = FALSE,
  theta = 0,
  phi = 30
)

fig <- plot_ly(
  x = output_vals[, 1],
  y = output_vals[, 3],
  z = output_vals[, 2],
  color = output_vals[, 4],
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 1.5)
) %>%
  layout(
    scene = list(
      camera = list(
        eye = list(
          x = 0,
          y = 2.25,
          z = 1
        )
      ),
      xaxis = list(title = "X<sub>1</sub>"),
      yaxis = list(title = "R"),
      zaxis = list(title = "X<sub>2</sub>")
    ),
    annotations = list(
      x = 1.075,
      y = 1.02,
      text = "Time",
      xref = "paper",
      yref = "paper",
      showarrow = FALSE
    )
  )
save_image(fig, "inconsistent_parameterization.png")
