library(plot3D)
library(pracma)
library(scatterplot3d)
library(tidyverse)

set.seed(1202)
theta <- runif(1000, min = 0, max = 2 * pi)
x <- cos(theta)
y <- sin(theta)

unit_df <- data.frame(theta, x, y) %>%
  arrange(theta)
angle <- cart2pol(as.matrix(unit_df[, 2:3], ncol = 2))


png("paper/figures/unit_circle_augmentation.png")
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
  zlab = "theta"
)
dev.off()

png("paper/figures/unit_circle_D2.png")
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

png("paper/figures/unit_circle_D3.png")
lines3D(
  x = unit_df$x,
  y = unit_df$y,
  z = angle[, 1],
  main = "Unit Circle, D = 3",
  xlab = "x",
  ylab = "y",
  zlab = "theta"
)
dev.off()

sim_example_case1 <- readRDS("simulations/case1/duration_01_interval_010_ampnoise_050_pernoise_050_n_1000_constant000_run_01.rds")
sim_example_case5 <- readRDS("simulations/case5/duration_01_interval_010_ampnoise_050_pernoise_000_n_10000_constant000_run_01.rds")
sim_example_case7 <- readRDS("simulations/case7/duration_01_interval_010_ampnoise_050_pernoise_050_n_1000_constant000_run_01.rds")

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

png("paper/figures/sim_case1.png")
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
    col = alpha("black", 0.4),
    bg = alpha("black", 0.4),
    xlim = c(-4, 4),
    ylim = c(-3, 3),
    pch = 21,
    main = paste0("Time = ", time_vals[time_idx])
  )
  points(
    x = temp_true$x,
    y = temp_true$y,
    col = "red",
    bg = "red",
    pch = 21
  )
  points(
    x = temp_lpme$x,
    y = temp_lpme$y,
    col = "blue",
    bg = "blue",
    pch = 21
  )
  points(
    x = temp_pme$x,
    y = temp_pme$y,
    col = "green",
    bg = "green",
    pch = 21
  )
  points(
    x = temp_princurve$x,
    y = temp_princurve$y,
    col = "purple",
    bg = "purple",
    pch = 21
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

png("paper/figures/sim_case5.png")
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
    xlim = c(-2, 2),
    ylim = c(0, 10),
    zlim = c(0, 20),
    pch = 21,
    col = alpha("black", 0.01),
    bg = alpha("black", 0.01),
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
    add = TRUE
  )
  scatter3D(
    x = temp_lpme$x,
    y = temp_lpme$y,
    z = temp_lpme$z,
    col = "blue",
    bg = "blue",
    pch = 20,
    add = TRUE
  )
  scatter3D(
    x = temp_pme$x,
    y = temp_pme$y,
    z = temp_pme$z,
    col = "green",
    bg = "green",
    pch = 20,
    add = TRUE
  )
  scatter3D(
    x = temp_princurve$x,
    y = temp_princurve$y,
    z = temp_princurve$z,
    col = "purple",
    bg = "purple",
    pch = 20,
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

png("paper/figures/sim_case7.png")
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
    pch = 21,
    col = alpha("black", 0.05),
    bg = alpha("black", 0.05),
    main = paste0("Time = ", time_vals[time_idx])
  )
  scatter3D(
    x = temp_true$x,
    y = temp_true$y,
    z = temp_true$z,
    col = alpha("red", 0.5),
    bg = alpha("red", 0.5),
    pch = 21,
    add = TRUE
  )
  scatter3D(
    x = temp_lpme$x,
    y = temp_lpme$y,
    z = temp_lpme$z,
    col = "blue",
    bg = "blue",
    pch = 21,
    add = TRUE
  )
  scatter3D(
    x = temp_pme$x,
    y = temp_pme$y,
    z = temp_pme$z,
    col = "green",
    bg = "green",
    pch = 21,
    add = TRUE
  )
  scatter3D(
    x = temp_princurve$x,
    y = temp_princurve$y,
    z = temp_princurve$z,
    col = "purple",
    bg = "purple",
    pch = 21,
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
