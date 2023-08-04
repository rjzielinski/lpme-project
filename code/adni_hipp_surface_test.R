library(cowplot)
library(gridGraphics)
library(pracma)
library(lubridate)
library(scatterplot3d)
library(tidyverse)
library(pme)
library(Rfast)

source("code/functions/calc_pme_est.R")
source("code/functions/calc_lpme_est.R")

create_scatter_mat <- function(mat, nrow = 1) {
  time_points <- unique(mat[, 1])
  xmin <- min(mat[, 2])
  xmax <- max(mat[, 2])
  ymin <- min(mat[, 3])
  ymax <- max(mat[, 3])
  zmin <- min(mat[, 4])
  zmax <- max(mat[, 4])

  # graphics::layout(matrix(c(1:length(time_points)), 1))

  plots <- list()
  for (time_idx in 1:length(time_points)) {
    time <- time_points[time_idx]
    scatterplot3d(
      x = mat[mat[, 1] == time, 2],
      y = mat[mat[, 1] == time, 3],
      z = mat[mat[, 1] == time, 4],
      main = paste0("Time = ", as.character(round(time, 3))),
      highlight.3d = TRUE,
      pch = 20,
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      zlim = c(zmin, zmax),
      xlab = "x",
      ylab = "y",
      zlab = "z",
      box = FALSE
    )
    plots[[time_idx]] <- recordPlot()
  }
  scatter_grid <- plot_grid(plotlist = plots, nrow = nrow)
  return(scatter_grid)
}

create_cross_section_matrix <- function(mats, slice_idx, slice_widths, nrow = 1) {
  time_vals <- unique(mats[[1]][, 1])
  slices <- list()
  for (time_idx in 1:length(time_vals)) {
    slice_list <- list()
    for (idx in 1:length(mats)) {
      mat <- mats[[idx]]
      temp_mat <- mat[mat[, 1] == time_vals[time_idx], ]
      n_pixels <- length(unique(temp_mat[, slice_idx]))
      slice_mean <- mean(temp_mat[, slice_idx])
      slice_list[[idx]] <- temp_mat[(temp_mat[, slice_idx] >= slice_mean - slice_widths[idx]) &
                                       (temp_mat[, slice_idx] <= slice_mean + slice_widths[idx]), ] %>%
        cbind(idx)
    }
    slices[[time_idx]] <- reduce(slice_list, rbind)
  }
  mat_sliced <- reduce(slices, rbind)
  # slice_pol <- cart2pol(mat_sliced[, 3:4])
  # mat_sliced <- cbind(mat_sliced[, 1:4], slice_pol)

  slice_df <- data.frame(mat_sliced)
  names(slice_df) <- c("time", "x", "y", "z", "theta", "phi", "r", "Source")
  slice_df <- slice_df %>%
    mutate(
      time = time * max_time,
      time = round(time, 2),
      time = as.factor(time),
      Source = ifelse(
        Source == 1,
        "Data",
        ifelse(
          Source == 2,
          "LPME",
          "PME"
        )
      ),
      opacity = case_when(
        Source == "Data" ~ 1,
        Source == "LPME" ~ 1,
        Source == "PME" ~ 1
      )
    )

  slice_df <- arrange(slice_df, Source)
  slice_df <- slice_df %>%
    mutate(
      Source = reorder(
        Source,
        X = c(
          rep(1, nrow(mat_sliced[mat_sliced[, 8] == 1, ])),
          rep(3, nrow(mat_sliced[mat_sliced[, 8] == 2, ])),
          rep(2, nrow(mat_sliced[mat_sliced[, 8] == 3, ]))
        )
      )
    )

  # slice_df <- slice_df %>%
  #   group_by(time) %>%
  #   arrange(time, phi)

  indices <- c("x", "y", "z")[2:4 != slice_idx]
  plt <- ggplot(arrange(slice_df, desc(as.character(Source)))) +
    geom_point(
      aes(
        x = .data[[indices[1]]],
        y = .data[[indices[2]]],
        color = Source
      ),
      size = 0.5
    ) +
    facet_wrap(vars(time), nrow = nrow)

  plt
}

lhipp_surface <- read_csv("data/adni_fsl_lhipp_surface_red.csv")

lhipp_surface <- lhipp_surface %>%
  mutate(
    scan_date = gsub("\\.0", "", scan_date),
    scan_date = ymd_hms(scan_date),
    scan_date = decimal_date(scan_date)
  )

lhipp_surface_centers <- lhipp_surface %>%
  group_by(
    patno,
    scan_date
  ) %>%
  summarize(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x)),
    max_y = max(abs(y)),
    max_z = max(abs(z))
  )

lhipp_bl <- lhipp_surface %>%
  group_by(patno) %>%
  arrange(scan_date) %>%
  summarize(
    time_bl = first(scan_date),
    max_time = max(scan_date)
  ) %>%
  mutate(duration = max_time - time_bl)

overlap <- 0.05

lhipp_surface <- lhipp_surface %>%
  full_join(lhipp_surface_centers, by = c("patno", "scan_date")) %>%
  full_join(lhipp_bl, by = "patno") %>%
  mutate(
    x = x - mean_x,
    y = y - mean_y,
    z = z - mean_z,
    range_x = max(x) - min(x),
    range_y = max(y) - min(y),
    range_z = max(z) - min(z)
  ) %>%
  mutate(
    partition1 = ifelse((x > (0 - (overlap * range_x))) & (y > (0 - (overlap * range_y))), TRUE, FALSE),
    partition2 = ifelse((x > (0 - (overlap * range_x))) & (y <= (0 + (overlap * range_y))), TRUE, FALSE),
    partition3 = ifelse((x <= (0 + (overlap * range_x))) & (y > (0 - (overlap * range_y))), TRUE, FALSE),
    partition4 = ifelse((x <= (0 + (overlap * range_x))) & (y <= (0 + (overlap * range_y))), TRUE, FALSE),
    partition_sum = partition1 + partition2 + partition3 + partition4,
    partition = case_when(
      partition1 == TRUE & partition_sum == 1 ~ 1,
      partition2 == TRUE & partition_sum == 1 ~ 2,
      partition3 == TRUE & partition_sum == 1 ~ 3,
      partition4 == TRUE & partition_sum == 1 ~ 4,
      partition_sum > 1 ~ 5
    ),
    time_from_bl = scan_date - time_bl
    # r = sqrt(x^2 + y^2 + z^2),
    # theta = atan(y / x),
    # phi = acos(z / r)
  ) %>%
  mutate(
    x = x / max_x,
    y = y / max_y,
    z = z / max_z,
    time_from_bl = time_from_bl / duration
  )

lhipp_surface_spherical <- lhipp_surface %>%
  dplyr::select(x, y, z) %>%
  as.matrix() %>%
  cart2sph() %>%
  as_tibble()

lhipp_surface <- bind_cols(
  lhipp_surface,
  lhipp_surface_spherical
)

lhipp_test <- lhipp_surface %>%
  filter(patno == "002_S_0413")

lhipp_test2 <- lhipp_surface %>%
  filter(patno == "002_S_0619")

lhipp_test3 <- lhipp_surface %>%
  filter(patno == "002_S_0685")

lhipp_test4 <- lhipp_surface %>%
  filter(patno == "002_S_0782")

# lhipp_plot1 <- plot_ly(
#   filter(lhipp_test4, time_from_bl == 0),
#   x = ~x,
#   y = ~y,
#   z = ~z,
#   type = "scatter3d",
#   mode = "markers",
#   marker = list(size = 3)
# )



lhipp_test_pt1 <- lhipp_test %>%
  filter(partition1 == TRUE)

lhipp_test_pt2 <- lhipp_test %>%
  filter(partition2 == TRUE)

lhipp_test_pt3 <- lhipp_test %>%
  filter(partition3 == TRUE)

lhipp_test_pt4 <- lhipp_test %>%
  filter(partition4 == TRUE)

lhipp_test_mat <- lhipp_test3 %>%
  dplyr::select(
    time_from_bl,
    x,
    y,
    z,
    theta,
    phi,
    r
  ) %>%
  as.matrix()

lhipp_test_bl <- lhipp_test_mat[lhipp_test_mat[, 1] == 0, -1]

lhipp_pt1_mat <- lhipp_test_pt1 %>%
  dplyr::select(
    time_from_bl,
    x,
    y,
    z
  ) %>%
  as.matrix()

lhipp_pt2_mat <- lhipp_test_pt2 %>%
  dplyr::select(
    time_from_bl,
    x,
    y,
    z
  ) %>%
  as.matrix()

lhipp_pt3_mat <- lhipp_test_pt3 %>%
  dplyr::select(
    time_from_bl,
    x,
    y,
    z
  ) %>%
  as.matrix()

lhipp_pt4_mat <- lhipp_test_pt4 %>%
  dplyr::select(
    time_from_bl,
    x,
    y,
    z
  ) %>%
  as.matrix()

# test_pme_pt1 <- pme(lhipp_pt1_mat[lhipp_pt1_mat[, 1] == 0, -1], 2)

# lpme_test_pt1 <- lpme(lhipp_pt1_mat, 2)
# lpme_test_pt2 <- lpme(lhipp_pt2_mat, 2)
# lpme_test_pt3 <- lpme(lhipp_pt3_mat, 2)
# lpme_test_pt4 <- lpme(lhipp_pt4_mat, 2)

# test_pme <- pme(lhipp_test_bl, d = 2, print_plots = TRUE)

time_vals <- unique(lhipp_test_mat[, 1])
# test_lpme <- lpme(lhipp_test_mat, d = 2, tuning.para.seq = exp(seq(-20, 5, 0.25)))
lhipp_test_mat5 <- lhipp_test_mat
lhipp_test_mat5[, 6] <- lhipp_test_mat5[, 6] * 2
test_lpme5 <- lpme(lhipp_test_mat5, d = 2, print_plots = TRUE)
lpme_vals <- calc_lpme_est(test_lpme, lhipp_test_mat)
# lpme_vals[, 1] <- sim_df[, 1]
pme_result <- list()
pme_vals <- list()
for (t in 1:length(time_vals)) {
  temp_data <- lhipp_test_mat5[lhipp_test_mat5[, 1] == time_vals[t], -1]
  pme_result[[t]] <- pme(temp_data, d = 2, verbose = "MSD")
  pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
}
pme_vals <- reduce(pme_vals, rbind)

plot_ly(
  x = lpme_preds5[, 2],
  y = lpme_preds5[, 3],
  z = lpme_preds5[, 4],
  frame = lpme_preds[, 1],
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3),
  name = "LPME"
) %>%
  add_markers(
    x = pme_preds[, 2],
    y = pme_preds[, 3],
    z = pme_preds[, 4],
    frame = pme_preds[, 1],
    name = "PME"
  ) %>%
  add_markers(
    x = lhipp_test_mat[, 2],
    y = lhipp_test_mat[, 3],
    z = lhipp_test_mat[, 4],
    frame = lhipp_test_mat[, 1],
    name = "Data"
  )

lpme_list <- list(
  lpme_result = test_lpme,
  df = lhipp_test_mat,
  lpme_vals = lpme_vals,
  pme_vals = pme_vals
)
saveRDS(lpme_list, "lhipp_test_results.RDS")

lpme_params <- test_lpme5$optimal_parameterization
lpme_param_extremes <- colMinsMaxs(lpme_params)
lpme_params_dim1 <- seq(
  lpme_param_extremes[1, 2],
  lpme_param_extremes[2, 2],
  length.out = 100
)
lpme_params_dim2 <- seq(
  lpme_param_extremes[1, 3],
  lpme_param_extremes[2, 3],
  length.out = 100
)
param_grid <- expand_grid(time_vals, lpme_params_dim1, lpme_params_dim2) %>%
  as.matrix()
lpme_preds5 <- map(
  1:nrow(param_grid),
  ~ test_lpme5$embedding_map(param_grid[.x, ])
) %>%
  reduce(rbind)

pme_params <- list()
pme_preds <- list()
for (time_idx in 1:length(time_vals)) {
  opt_params <- pme_result[[time_idx]]$parameterization[[which.min(pme_result[[time_idx]]$MSD)]]
  param_extremes <- colMinsMaxs(opt_params)
  param_dim1 <- seq(
    param_extremes[1, 1],
    param_extremes[2, 1],
    length.out = 100
  )
  param_dim2 <- seq(
    param_extremes[1, 2],
    param_extremes[2, 2],
    length.out = 100
  )
  pme_params[[time_idx]] <- cbind(
    time_vals[time_idx],
    as.matrix(expand_grid(param_dim1, param_dim2))
  )
  pme_preds[[time_idx]] <- map(
    1:nrow(pme_params[[time_idx]]),
    ~ pme_result[[time_idx]]$embedding_map(pme_params[[time_idx]][.x, -1])
  ) %>%
    reduce(rbind)
  pme_preds[[time_idx]] <- cbind(time_vals[time_idx], pme_preds[[time_idx]])
}
pme_preds <- reduce(pme_preds, rbind)

create_cross_section_matrix(lhipp_test_mat, 2, slice_width = 0.005, nrow = 2)
create_cross_section_matrix(lpme_preds, 2, slice_width = 0.003, nrow = 2)
create_cross_section_matrix(pme_preds, 2, slice_width = 0.003, nrow = 2)

### Next step: glue estimated manifolds together to estimate full surface

max_time <- lhipp_surface %>%
  filter(patno == "002_S_0685") %>%
  select(duration) %>%
  unique() %>%
  as.numeric()

time_vals <- unique(lhipp_test_mat[, 1])
png("paper/figures/adni_lhipp.png", res = 500, height = 2000, width = 3000)
par(oma = c(4, 1, 1, 1), mfrow = c(2, 3), mar = c(2, 2, 1, 1))
for (time_idx in 1:length(time_vals)) {
  temp_data <- lhipp_test_mat[lhipp_test_mat[, 1] == time_vals[time_idx], ]
  temp_lpme <- lpme_preds[lpme_preds[, 1] == time_vals[time_idx], ]
  temp_pme <- pme_preds[pme_vals[, 1] == time_vals[time_idx], ]

  plt <- scatter3D(
    x = temp_data[, 2],
    y = temp_data[, 3],
    z = temp_data[, 4],
    xlab = "x",
    ylab = "y",
    zlab = "z",
    xlim = c(-0.2, 0.2),
    ylim = c(-0.2, 0.2),
    zlim = c(-0.2, 0.2),
    pch = 21,
    col = alpha("black", 0.1),
    bg = alpha("black", 0.1),
    cex = 0.33,
    main = paste0("Time = ", round(time_vals[time_idx] * max_time, 2)),
    theta = 225,
    phi = 60
  )
  scatter3D(
    x = temp_lpme[, 2],
    y = temp_lpme[, 3],
    z = temp_lpme[, 4],
    col = alpha("blue", 0.1),
    bg = alpha("blue", 0.1),
    pch = 21,
    cex = 0.33,
    add = TRUE
  )
  scatter3D(
    x = temp_pme[, 2],
    y = temp_pme[, 3],
    z = temp_pme[, 4],
    col = alpha("green", 0.1),
    bg = alpha("green", 0.1),
    pch = 21,
    cex = 0.33,
    add = TRUE
  )
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "l", bty = "n", xaxt = "n", yaxt = "n")
legend(
  "bottom",
  c("Data", "LPME", "PME"),
  col = c("black", "blue", "green"),
  # fill = c("black", "red", "blue", "green", "purple"),
  lwd = 5,
  xpd = TRUE,
  horiz = TRUE,
  cex = 1,
  bty = "n"
)
dev.off()

lhipp_maxes <- lhipp_surface %>%
  filter(patno == "002_S_0685") %>%
  group_by(time_from_bl) %>%
  summarize(
    max_x = unique(max_x),
    max_y = unique(max_y),
    max_z = unique(max_z)
  )

lhipp_time_idx <- as.factor(lpme_vals[, 1]) %>%
  as.numeric()

lpme_vals_mod <- lpme_vals
lpme_vals_mod[, 2] <- round(lpme_vals_mod[, 2] * lhipp_maxes[lhipp_time_idx, 2]) %>%
  as.matrix()
lpme_vals_mod[, 3] <- round(lpme_vals_mod[, 3] * lhipp_maxes[lhipp_time_idx, 3]) %>%
  as.matrix()
lpme_vals_mod[, 4] <- round(lpme_vals_mod[, 4] * lhipp_maxes[lhipp_time_idx, 4]) %>%
  as.matrix()

