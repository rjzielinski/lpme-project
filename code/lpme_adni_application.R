library(lubridate)
library(plotly)
library(plot3D)
library(pme)
library(pracma)
library(RColorBrewer)
library(Rfast)
library(tidyverse)

source("code/functions/calc_pme_est.R")
source("code/functions/calc_lpme_est_params.R")
source("code/functions/estimate_volume.R")
source("code/functions/interior_identification.R")

# function for creating cross-section plots
create_cross_section_matrix <- function(mats, slice_idx, slice_widths, nrow = 1, lhipp) {
  time_vals <- unique(mats[[1]][, 1])
  slices <- list()
  for (time_idx in 1:length(time_vals)) {
    slice_list <- list()
    for (idx in 1:length(mats)) {
      mat <- mats[[idx]]
      temp_mat <- mat[mat[, 1] == time_vals[time_idx], ]
      n_pixels <- length(unique(temp_mat[, slice_idx]))
      slice_mean <- mean(temp_mat[, slice_idx])
      # in the PME- and LPME-generated value matrices, voxel values are not all
      # consistent as in the data
      # to approximate a cross-section through the middle of the structure, we
      # consider all voxel values within a certain range of the mean voxel value.
      slice_list[[idx]] <- temp_mat[(temp_mat[, slice_idx] >= slice_mean - slice_widths[idx]) &
                                      (temp_mat[, slice_idx] <= slice_mean + slice_widths[idx]), ] %>%
        cbind(idx)
    }
    slices[[time_idx]] <- reduce(slice_list, rbind)
  }
  mat_sliced <- reduce(slices, rbind)

  slice_df <- data.frame(mat_sliced)
  names(slice_df) <- c("time", "x", "y", "z", "theta", "phi", "r", "Source")
  slice_df <- slice_df %>%
    mutate(
      time = time * unique(lhipp$duration),
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
      # code to change opacity by data source if needed
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

# lhipp_surface <- read_csv("data/adni_fsl_lhipp_surface_red.csv")
lhipp_surface <- read_csv("~/Documents/brown/research/lpme-project/data/adni_fsl_lhipp_surface_red.csv")
lthal_surface <- read_csv("~/Documents/brown/research/lpme-project/data/adni_fsl_lthal_surface.csv")

# process scan dates
lhipp_surface <- lhipp_surface %>%
  mutate(
    scan_date = gsub("\\.0", "", scan_date),
    scan_date = ymd_hms(scan_date),
    scan_date = decimal_date(scan_date)
  )
lthal_surface <- lthal_surface %>%
  mutate(
    scan_date = gsub("\\.0", "", scan_date),
    scan_date = ymd_hms(scan_date),
    scan_date = decimal_date(scan_date)
  )

# lhipp_surface_spherical <- lhipp_surface %>%
#   dplyr::select(x, y, z) %>%
#   as.matrix() %>%
#   cart2sph() %>%
#   as_tibble()

# lhipp_surface <- bind_cols(lhipp_surface, lhipp_surface_spherical)

lhipp_surface_centers <- lhipp_surface %>%
  group_by(patno, scan_date) %>%
  summarize(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    # mean_theta = mean(theta),
    # mean_phi = mean(phi),
    # mean_r = mean(r),
    # max_x = max(abs(x - mean_x)),
    # max_y = max(abs(y - mean_y)),
    # max_z = max(abs(z - mean_z)),
    max_x = max(abs(x)),
    max_y = max(abs(y)),
    max_z = max(abs(z)),
    sd_x = sd(x),
    sd_y = sd(y),
    sd_z = sd(z)
    # max_theta = max(abs(theta - mean_theta)),
    # max_phi = max(abs(phi - mean_phi)),
    # max_r = max(abs(r - mean_r))
    # max_theta = max(abs(theta)),
    # max_phi = max(abs(phi)),
    # max_r = max(abs(r)),
    # sd_theta = sd(theta),
    # sd_phi = sd(phi),
    # sd_r = sd(r)
  )

lthal_surface_centers <- lthal_surface %>%
  group_by(patno, scan_date) %>%
  summarize(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    # mean_theta = mean(theta),
    # mean_phi = mean(phi),
    # mean_r = mean(r),
    # max_x = max(abs(x - mean_x)),
    # max_y = max(abs(y - mean_y)),
    # max_z = max(abs(z - mean_z)),
    max_x = max(abs(x)),
    max_y = max(abs(y)),
    max_z = max(abs(z)),
    sd_x = sd(x),
    sd_y = sd(y),
    sd_z = sd(z)
    # max_theta = max(abs(theta - mean_theta)),
    # max_phi = max(abs(phi - mean_phi)),
    # max_r = max(abs(r - mean_r))
    # max_theta = max(abs(theta)),
    # max_phi = max(abs(phi)),
    # max_r = max(abs(r)),
    # sd_theta = sd(theta),
    # sd_phi = sd(phi),
    # sd_r = sd(r)
  )

lhipp_bl <- lhipp_surface %>%
  group_by(patno) %>%
  arrange(scan_date) %>%
  summarize(
    time_bl = first(scan_date),
    max_time = max(scan_date)
  ) %>%
  mutate(duration = max_time - time_bl)

lthal_bl <- lthal_surface %>%
  group_by(patno) %>%
  arrange(scan_date) %>%
  summarize(
    time_bl = first(scan_date),
    max_time = max(scan_date)
  ) %>%
  mutate(duration = max_time - time_bl)

# center and standardize data
lhipp_surface <- lhipp_surface %>%
  full_join(lhipp_surface_centers, by = c("patno", "scan_date")) %>%
  full_join(lhipp_bl, by = "patno") %>%
  mutate(
    range_x = max(x) - min(x),
    range_y = max(y) - min(y),
    range_z = max(z) - min(z),
    x = (x - mean_x) / max_x,
    y = (y - mean_y) / max_y,
    z = (z - mean_z) / max_z,
    # theta = (theta - mean_theta) / max_theta,
    # phi = (phi - mean_phi) / max_phi,
    # r = (r - mean_r) / max_r,
    # x = (x - mean_x) / sd_x,
    # y = (y - mean_y) / sd_y,
    # z = (z - mean_z) / sd_z,
    # theta = (theta - mean_theta) / sd_theta,
    # phi = (phi - mean_phi) / sd_phi,
    # r = (r - mean_r) / sd_r,
    time_from_bl = (scan_date - time_bl) / duration
  )

lthal_surface <- lthal_surface %>%
  full_join(lthal_surface_centers, by = c("patno", "scan_date")) %>%
  full_join(lthal_bl, by = "patno") %>%
  mutate(
    range_x = max(x) - min(x),
    range_y = max(y) - min(y),
    range_z = max(z) - min(z),
    x = (x - mean_x) / max_x,
    y = (y - mean_y) / max_y,
    z = (z - mean_z) / max_z,
    # theta = (theta - mean_theta) / max_theta,
    # phi = (phi - mean_phi) / max_phi,
    # r = (r - mean_r) / max_r,
    # x = (x - mean_x) / sd_x,
    # y = (y - mean_y) / sd_y,
    # z = (z - mean_z) / sd_z,
    # theta = (theta - mean_theta) / sd_theta,
    # phi = (phi - mean_phi) / sd_phi,
    # r = (r - mean_r) / sd_r,
    time_from_bl = (scan_date - time_bl) / duration
  )

 lhipp_surface_spherical <- lhipp_surface %>%
   dplyr::select(x, y, z) %>%
   as.matrix() %>%
   cart2sph() %>%
   as_tibble()

 lthal_surface_spherical <- lthal_surface %>%
   dplyr::select(x, y, z) %>%
   as.matrix() %>%
   cart2sph() %>%
   as_tibble()

 lhipp_surface <- bind_cols(lhipp_surface, lhipp_surface_spherical)
 lthal_surface <- bind_cols(lthal_surface, lthal_surface_spherical)

 # lhipp_surface_spherical_centers <- lhipp_surface %>%
 #   group_by(patno, scan_date) %>%
 #   summarize(
 #     mean_theta = mean(theta),
 #     mean_phi = mean(phi),
 #     mean_r = mean(r),
 #     max_theta = max(abs(theta)),
 #     max_phi = max(abs(phi)),
 #     max_r = max(abs(r)),
 #     sd_theta = sd(theta),
 #     sd_phi = sd(phi),
 #     sd_r = sd(r)
 #   )
 #
 # lhipp_surface <- full_join(
 #   lhipp_surface,
 #   lhipp_surface_spherical_centers,
 #   by = c("patno", "scan_date")
 # ) %>%
 #   mutate(
 #   # theta = (theta - mean_theta) / max_theta,
 #   # phi = (phi - mean_phi) / max_phi,
 #   # r = (r - mean_r) / max_r
 #   theta = (theta - mean_theta) / sd_theta,
 #   phi = (phi - mean_phi) / sd_phi,
 #   r = (r - mean_r) / sd_r
 #   )

lhipp <- lhipp_surface %>%
  filter(patno == "002_S_0685")

lthal <- lthal_surface %>%
  filter(patno == "002_S_0685")

lhipp_mat <- lhipp %>%
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

lthal_mat <- lthal %>%
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

time_vals <- unique(lhipp_mat[, 1])

# set.seed(1036)
set.seed(10283)
lhipp_lpme_result_isomap <- lpme(
  lhipp_mat,
  d = 2,
  gamma = exp(-15:5),
  lambda = exp(-12:5),
  initialization_algorithm = "isomap",
  init_type = "subsampling",
  verbose = TRUE,
  print_plots = TRUE
)
lhipp_lpme_isomap_est <- calc_lpme_est_params(lhipp_lpme_result_isomap, lhipp_mat)
lhipp_lpme_isomap_vals <- lhipp_lpme_isomap_est$results
lhipp_lpme_isomap_params <- lhipp_lpme_isomap_est$params

# lhipp_lpme_result_dm <- lpme(
#   lhipp_mat,
#   d = 2,
#   gamma = exp(-15:5),
#   lambda = exp(-15:5),
#   initialization_algorithm = "diffusion_maps",
#   init_type = "subsampling",
#   verbose = TRUE,
#   print_plots = TRUE
# )
# lhipp_lpme_dm_est <- calc_lpme_est_params(lhipp_lpme_result_dm, lhipp_mat)
# lhipp_lpme_dm_vals <- lhipp_lpme_dm_est$results
# lhipp_lpme_dm_params <- lhipp_lpme_dm_est$params
#
# lhipp_lpme_result_laplacian <- lpme(
#   lhipp_mat,
#   d = 2,
#   gamma = exp(-15:5),
#   lambda = exp(-15:5),
#   initialization_algorithm = "laplacian_eigenmaps",
#   init_type = "subsampling",
#   verbose = TRUE,
#   print_plots = TRUE
# )
# lhipp_lpme_laplacian_est <- calc_lpme_est_params(lhipp_lpme_result_laplacian, lhipp_mat)
# lhipp_lpme_laplacian_vals <- lhipp_lpme_laplacian_est$results
# lhipp_lpme_laplacian_params <- lhipp_lpme_laplacian_est$params

lhipp_pme_result <- list()
lhipp_pme_vals <- list()
for (t in 1:length(time_vals)) {
  temp_data <- lhipp_mat[lhipp_mat[, 1] == time_vals[t], -1]
  lhipp_pme_result[[t]] <- pme(temp_data, d = 2, lambda = exp(-12:5), verbose = TRUE)
  lhipp_pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(lhipp_pme_result[[t]], temp_data))
}
lhipp_pme_vals <- reduce(lhipp_pme_vals, rbind)

write.csv(lhipp_mat, "results/adni_lhipp_mat.csv")
write.csv(lhipp_lpme_isomap_vals, "results/adni_lhipp_lpme_isomap_vals.csv")
write.csv(lhipp_lpme_dm_vals, "results/adni_lhipp_lpme_dm_vals.csv")
write.csv(lhipp_lpme_laplacian_vals, "results/adni_lhipp_lpme_laplacian_vals.csv")
write.csv(lhipp_pme_vals, "results/adni_lhipp_pme_vals.csv")

p <- plot_ly(
  x = lhipp_mat[, 2],
  y = lhipp_mat[, 3],
  z = lhipp_mat[, 4],
  frame = lhipp_mat[, 1],
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3),
  name = "Data"
) %>%
  add_markers(
    x = lhipp_lpme_isomap_vals[, 2],
    y = lhipp_lpme_isomap_vals[, 3],
    z = lhipp_lpme_isomap_vals[, 4],
    frame = lhipp_lpme_isomap_vals[, 1],
    name = "LPME"
  ) %>%
  add_markers(
    x = lhipp_pme_vals[, 2],
    y = lhipp_pme_vals[, 3],
    z = lhipp_pme_vals[, 4],
    frame = lhipp_pme_vals[, 1],
    name = "PME"
  )

lhipp_lpme_plots <- list()
lhipp_pme_plots <- list()
lhipp_data_plots <- list()

colors <- brewer.pal(3, "Set1")

png(
  "results/adni_plots/adni_lhipp_data_plot.png",
  res = 1000,
  height = 15000,
  width = 3500
)
par(oma = c(4, 1, 1, 1), mfrow = c(6, 1), mar = c(2, 2, 1, 1))

for (t in 1:length(time_vals)) {
  temp_data <- lhipp_mat[lhipp_mat[, 1] == time_vals[t], ]
  scatter3D(
    x = temp_data[, 2],
    y = temp_data[, 3],
    z = temp_data[, 4],
    pch = 20,
    col = alpha(colors[1], 0.5),
    theta = 220,
    ticktype = "detailed",
    main = paste0("Time = ", round((time_vals * max(lhipp$duration))[t], 2)),
    xlim = c(-0.15, 0.15),
    ylim = c(-0.15, 0.15),
    zlim = c(-0.15, 0.15)
  )
}
dev.off()

png(
  "results/adni_plots/adni_lhipp_lpme_isomap_plot.png",
  res = 1000,
  height = 15000,
  width = 3500
)
par(oma = c(4, 1, 1, 1), mfrow = c(6, 1), mar = c(2, 2, 1, 1))

for (t in 1:length(time_vals)) {
  temp_lpme <- lhipp_lpme_isomap_vals[lhipp_lpme_isomap_vals[, 1] == time_vals[t], ]
  scatter3D(
    x = temp_lpme[, 2],
    y = temp_lpme[, 3],
    z = temp_lpme[, 4],
    pch = 20,
    col = alpha(colors[2], 0.5),
    theta = 220,
    ticktype = "detailed",
    main = paste0("Time = ", round((time_vals * max(lhipp$duration))[t], 2)),
    xlim = c(-0.15, 0.15),
    ylim = c(-0.15, 0.15),
    zlim = c(-0.15, 0.15)
  )
}
dev.off()

png(
  "results/adni_plots/adni_lhipp_lpme_dm_plot.png",
  res = 1000,
  height = 15000,
  width = 3500
)
par(oma = c(4, 1, 1, 1), mfrow = c(6, 1), mar = c(2, 2, 1, 1))

for (t in 1:length(time_vals)) {
  temp_lpme <- lhipp_lpme_dm_vals[lhipp_lpme_dm_vals[, 1] == time_vals[t], ]
  scatter3D(
    x = temp_lpme[, 2],
    y = temp_lpme[, 3],
    z = temp_lpme[, 4],
    pch = 20,
    col = alpha(colors[2], 0.5),
    theta = 220,
    ticktype = "detailed",
    main = paste0("Time = ", round((time_vals * max(lhipp$duration))[t], 2)),
    xlim = c(-0.15, 0.15),
    ylim = c(-0.15, 0.15),
    zlim = c(-0.15, 0.15)
  )
}
dev.off()

png(
  "results/adni_plots/adni_lhipp_lpme_laplacian_plot.png",
  res = 1000,
  height = 15000,
  width = 3500
)
par(oma = c(4, 1, 1, 1), mfrow = c(6, 1), mar = c(2, 2, 1, 1))

for (t in 1:length(time_vals)) {
  temp_lpme <- lhipp_lpme_laplacian_vals[lhipp_lpme_laplacian_vals[, 1] == time_vals[t], ]
  scatter3D(
    x = temp_lpme[, 2],
    y = temp_lpme[, 3],
    z = temp_lpme[, 4],
    pch = 20,
    col = alpha(colors[2], 0.5),
    theta = 220,
    ticktype = "detailed",
    main = paste0("Time = ", round((time_vals * max(lhipp$duration))[t], 2)),
    xlim = c(-0.15, 0.15),
    ylim = c(-0.15, 0.15),
    zlim = c(-0.15, 0.15)
  )
}
dev.off()

png(
  "results/adni_plots/adni_lhipp_pme_plot.png",
  res = 1000,
  height = 15000,
  width = 3500
)
par(oma = c(4, 1, 1, 1), mfrow = c(6, 1), mar = c(2, 2, 1, 1))

for (t in 1:length(time_vals)) {
  temp_pme <- lhipp_pme_vals[lhipp_pme_vals[, 1] == time_vals[t], ]
  scatter3D(
    x = temp_pme[, 2],
    y = temp_pme[, 3],
    z = temp_pme[, 4],
    pch = 20,
    col = alpha(colors[3], 0.5),
    theta = 220,
    ticktype = "detailed",
    main = paste0("Time = ", round((time_vals * max(lhipp$duration))[t], 2)),
    xlim = c(-0.15, 0.15),
    ylim = c(-0.15, 0.15),
    zlim = c(-0.15, 0.15)
  )
}
dev.off()

lhipp_cross_section <- create_cross_section_matrix(
  list(lhipp_mat, lhipp_lpme_isomap_vals, lhipp_pme_vals),
  2,
  c(0.005, 0.005, 0.005),
  nrow = 2,
  lhipp = lhipp
)
ggsave(
  "results/adni_plots/adni_lhipp_cross_section.png",
  plot = lhipp_cross_section,
  height = 5000,
  width = 7000,
  units = "px",
  dpi = 1250
)

#### VOLUME ESTIMATION

# Voxel dimension: 1.2mm x 0.9375mm x 0.9375mm
# Voxel volume: 1.054688 mm^3

lhipp_rescaled <- lhipp %>%
  mutate(
    x_rescaled = round((x * max_x)),
    y_rescaled = round((y * max_y)),
    z_rescaled = round((z * max_z))
  )
candidate_x <- seq(min(lhipp_rescaled$x_rescaled), max(lhipp_rescaled$x_rescaled), 1)
candidate_y <- seq(min(lhipp_rescaled$y_rescaled), max(lhipp_rescaled$y_rescaled), 1)
candidate_z <- seq(min(lhipp_rescaled$z_rescaled), max(lhipp_rescaled$z_rescaled), 1)

candidate_voxels <- expand_grid(candidate_x, candidate_y, candidate_z)

lhipp_pme_volumes <- vector(mode = "numeric", length = length(time_vals))
lhipp_lpme_volumes <- vector(mode = "numeric", length = length(time_vals))
lhipp_pme_interior_plots <- list()
lhipp_lpme_interior_plots <- list()

for (time_idx in 1:length(time_vals)) {
  temp_max_x <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$max_x)
  temp_max_y <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$max_y)
  temp_max_z <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$max_z)

  temp_candidates <- candidate_voxels %>%
    mutate(
      x_scaled = candidate_x / temp_max_x,
      y_scaled = candidate_y / temp_max_y,
      z_scaled = candidate_z / temp_max_z
    )

  temp_pme <- lhipp_pme_result[[time_idx]]

  temp_pme_embedding <- temp_pme$embedding_map
  temp_pme_coefs <- temp_pme$coefs[[which.min(temp_pme$MSD)]]
  temp_pme_params <- temp_pme$parameterization[[which.min(temp_pme$MSD)]]

  temp_lpme_coefs <- lhipp_lpme_result_isomap$sol_coef_functions[[which.min(lhipp_lpme_result_isomap$msd)]](time_vals[time_idx])
  temp_lpme_params <- lhipp_lpme_result_isomap$parameterization_list[[which.min(lhipp_lpme_result_isomap$msd)]]
  temp_lpme_params <- temp_lpme_params[temp_lpme_params[, 1] == time_vals[time_idx], -1]
  lpme_n_knots <- nrow(temp_lpme_params)
  d <- ncol(temp_lpme_params)
  coef_mat <- matrix(
    temp_lpme_coefs,
    lpme_n_knots + d + 1,
    byrow = TRUE
  )

  temp_lpme_embedding <- function(r) {
    t(coef_mat[1:lpme_n_knots, ]) %*% pme::etaFunc(r, temp_lpme_params, 4 - d) +
      t(coef_mat[(lpme_n_knots + 1):(lpme_n_knots + d + 1), ]) %*% matrix(c(1, r), ncol = 1)
  }

  lhipp_interior_voxel_pme <- vector(length = nrow(temp_candidates))
  lhipp_interior_voxel_lpme <- vector(length = nrow(temp_candidates))
  for (row_idx in 1:nrow(temp_candidates)) {
    lhipp_interior_voxel_pme[row_idx] <- interior_identification(
      temp_pme_embedding,
      temp_pme_coefs,
      temp_pme_params,
      unlist(temp_candidates[row_idx, 4:6]),
      c(0, 0, 0)
    )
    lhipp_interior_voxel_lpme[row_idx] <- interior_identification(
      temp_lpme_embedding,
      coef_mat,
      temp_lpme_params,
      unlist(temp_candidates[row_idx, 4:6]),
      c(0, 0, 0)
    )
  }

  lhipp_interior_points_pme <- temp_candidates[lhipp_interior_voxel_pme, ]
  lhipp_interior_points_lpme <- temp_candidates[lhipp_interior_voxel_lpme, ]
  p <- plot_ly(
    lhipp_interior_points_pme,
    x = ~candidate_x,
    y = ~candidate_y,
    z = ~candidate_z,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  )

  p_lpme <- plot_ly(
    lhipp_interior_points_lpme,
    x = ~candidate_x,
    y = ~candidate_y,
    z = ~candidate_z,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  )

  lhipp_pme_interior_plots[[time_idx]] <- p
  lhipp_lpme_interior_plots[[time_idx]] <- p_lpme
  lhipp_pme_volumes[time_idx] <- 1.054688 * nrow(lhipp_interior_points_pme)
  lhipp_lpme_volumes[time_idx] <- 1.054688 * nrow(lhipp_interior_points_lpme)
}


# ##### Interior Identification Test on Spherical Simulated Data
#
# # sim_df and pme_result are taken from sim_error_case9() in lpme_simulations.R
# test_pme <- pme_result[[1]]
# test_df <- sim_df[sim_df[, 1] == 0, -1]
#
# candidate_x <- seq(min(test_df[, 1]), max(test_df[, 1]), length.out = 25)
# candidate_y <- seq(min(test_df[, 2]), max(test_df[, 2]), length.out = 25)
# candidate_z <- seq(min(test_df[, 3]), max(test_df[, 3]), length.out = 25)
#
# candidate_voxels <- expand_grid(candidate_x, candidate_y, candidate_z)
#
# interior_voxels <- map(
#   1:nrow(candidate_voxels),
#   ~ interior_identification(
#     test_pme$embedding_map,
#     test_pme$coefs[[which.min(test_pme$MSD)]],
#     test_pme$parameterization[[which.min(test_pme$MSD)]],
#     unlist(candidate_voxels[.x, ]),
#     c(0, 0, 0.2)
#   )
# ) %>%
#   reduce(c)
#
# interior_points <- candidate_voxels[interior_voxels, ]
# plot_ly(
#   interior_points,
#   x = ~candidate_x,
#   y = ~candidate_y,
#   z = ~candidate_z,
#   type = "scatter3d",
#   mode = "markers",
#   marker = list(size = 3)
# )

# Begin by rescaling data and estimated surfaces to full size

# lhipp_lpme_vals <- lhipp_lpme_isomap_vals
#
# lhipp_data_rescaled <- list()
# lhipp_lpme_rescaled <- list()
# lhipp_pme_rescaled <- list()
#
# for (time_idx in 1:length(time_vals)) {
#   temp_data <- lhipp_mat[lhipp_mat[, 1] == time_vals[time_idx], ]
#   temp_lpme <- lhipp_lpme_vals[lhipp_lpme_vals[, 1] == time_vals[time_idx], ]
#   temp_pme <- lhipp_pme_vals[lhipp_pme_vals[, 1] == time_vals[time_idx], ]
#
#   temp_max_x <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$max_x)
#   temp_max_y <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$max_y)
#   temp_max_z <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$max_z)
#   # temp_sd_x <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$sd_x)
#   # temp_sd_y <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$sd_y)
#   # temp_sd_z <- unique(lhipp[lhipp$time_from_bl == time_vals[time_idx], ]$sd_z)
#
#   temp_data_rescaled <- temp_data[, -(5:7)]
#   temp_lpme_rescaled <- temp_lpme[, -(5:7)]
#   temp_pme_rescaled <- temp_pme[, -(5:7)]
#
#   temp_data_rescaled[, 2] <- round((temp_data_rescaled[, 2] * temp_max_x) / 1.2)
#   temp_data_rescaled[, 3] <- round((temp_data_rescaled[, 3] * temp_max_y) / 0.9375)
#   temp_data_rescaled[, 4] <- round((temp_data_rescaled[, 4] * temp_max_z) / 0.9375)
#
#   temp_lpme_rescaled[, 2] <- round((temp_lpme_rescaled[, 2] * temp_max_x) / 1.2)
#   temp_lpme_rescaled[, 3] <- round((temp_lpme_rescaled[, 3] * temp_max_y) / 0.9375)
#   temp_lpme_rescaled[, 4] <- round((temp_lpme_rescaled[, 4] * temp_max_z) / 0.9375)
#
#   temp_pme_rescaled[, 2] <- round((temp_pme_rescaled[, 2] * temp_max_x) / 1.2)
#   temp_pme_rescaled[, 3] <- round((temp_pme_rescaled[, 3] * temp_max_y) / 0.9375)
#   temp_pme_rescaled[, 4] <- round((temp_pme_rescaled[, 4] * temp_max_z) / 0.9375)
#
#   # temp_data_rescaled[, 2] <- round((temp_data_rescaled[, 2] * temp_sd_x) / 1.2)
#   # temp_data_rescaled[, 3] <- round((temp_data_rescaled[, 3] * temp_sd_y) / 0.9375)
#   # temp_data_rescaled[, 4] <- round((temp_data_rescaled[, 4] * temp_sd_z) / 0.9375)
#   #
#   # temp_lpme_rescaled[, 2] <- round((temp_lpme_rescaled[, 2] * temp_sd_x) / 1.2)
#   # temp_lpme_rescaled[, 3] <- round((temp_lpme_rescaled[, 3] * temp_sd_y) / 0.9375)
#   # temp_lpme_rescaled[, 4] <- round((temp_lpme_rescaled[, 4] * temp_sd_z) / 0.9375)
#   #
#   # temp_pme_rescaled[, 2] <- round((temp_pme_rescaled[, 2] * temp_sd_x) / 1.2)
#   # temp_pme_rescaled[, 3] <- round((temp_pme_rescaled[, 3] * temp_sd_y) / 0.9375)
#   # temp_pme_rescaled[, 4] <- round((temp_pme_rescaled[, 4] * temp_sd_z) / 0.9375)
#
#   lhipp_data_rescaled[[time_idx]] <- temp_data_rescaled
#   lhipp_lpme_rescaled[[time_idx]] <- temp_lpme_rescaled
#   lhipp_pme_rescaled[[time_idx]] <- temp_pme_rescaled
# }
#
# lhipp_data_rescaled_full <- reduce(lhipp_data_rescaled, rbind)
# lhipp_lpme_rescaled_full <- reduce(lhipp_lpme_rescaled, rbind)
# lhipp_pme_rescaled_full <- reduce(lhipp_pme_rescaled, rbind)
#
# voxel_vol <- 1.2 * 0.9375 * 0.9375
#
# lhipp_est_data_volume <- map(
#   lhipp_data_rescaled,
#   ~ estimate_volume(.x[, -1], voxel_vol)
# )
#
# lhipp_est_data_volume_points <- map(
#   lhipp_est_data_volume,
#   ~ .x[[2]]
# )
#
# lhipp_est_data_volume <- map(
#   lhipp_est_data_volume,
#   ~ .x[[1]]
# ) %>%
#   reduce(c)
#
# lhipp_est_lpme_volume <- map(
#   lhipp_lpme_rescaled,
#   ~ estimate_volume(.x[, -1], voxel_vol)
# )
#
# lhipp_est_lpme_volume_points <- map(
#   lhipp_est_lpme_volume,
#   ~ .x[[2]]
# )
#
# lhipp_est_lpme_volume <- map(
#   lhipp_est_lpme_volume,
#   ~ .x[[1]]
# ) %>%
#   reduce(c)
#
# lhipp_est_pme_volume <- map(
#   lhipp_pme_rescaled,
#   ~ estimate_volume(.x[, -1], voxel_vol)
# )
#
# lhipp_est_pme_volume_points <- map(
#   lhipp_est_pme_volume,
#   ~ .x[[2]]
# )
#
# lhipp_est_pme_volume <- map(
#   lhipp_est_pme_volume,
#   ~ .x[[1]]
# ) %>%
#   reduce(c)
#
# lhipp_est_data_volume_points_full <- map(
#   1:length(lhipp_est_data_volume_points),
#   ~ cbind(time_vals[.x], lhipp_est_data_volume_points[[.x]])
# ) %>%
#   reduce(rbind)
#
#
# lhipp_est_lpme_volume_points_full <- map(
#   1:length(lhipp_est_lpme_volume_points),
#   ~ cbind(time_vals[.x], lhipp_est_lpme_volume_points[[.x]])
# ) %>%
#   reduce(rbind)
#
# lhipp_est_pme_volume_points_full <- map(
#   1:length(lhipp_est_pme_volume_points),
#   ~ cbind(time_vals[.x], lhipp_est_pme_volume_points[[.x]])
# ) %>%
#   reduce(rbind)
#


set.seed(61829)
lthal_lpme_result_isomap <- lpme(
  lthal_mat,
  d = 2,
  gamma = exp(-15:5),
  lambda = exp(-12:5),
  initialization_algorithm = "isomap",
  init_type = "subsampling",
  verbose = TRUE,
  print_plots = TRUE
)
lthal_lpme_isomap_est <- calc_lpme_est_params(lthal_lpme_result_isomap, lthal_mat)
lthal_lpme_isomap_vals <- lthal_lpme_isomap_est$results
lthal_lpme_isomap_params <- lthal_lpme_isomap_est$params

# lthal_lpme_result_dm <- lpme(
#   lthal_mat,
#   d = 2,
#   gamma = exp(-15:5),
#   lambda = exp(-15:5),
#   initialization_algorithm = "diffusion_maps",
#   init_type = "subsampling",
#   verbose = TRUE,
#   print_plots = TRUE
# )
# lthal_lpme_dm_est <- calc_lpme_est_params(lthal_lpme_result_dm, lthal_mat)
# lthal_lpme_dm_vals <- lthal_lpme_dm_est$results
# lthal_lpme_dm_params <- lthal_lpme_dm_est$params
#
# lthal_lpme_result_laplacian <- lpme(
#   lthal_mat,
#   d = 2,
#   gamma = exp(-15:5),
#   lambda = exp(-15:5),
#   initialization_algorithm = "laplacian_eigenmaps",
#   init_type = "subsampling",
#   verbose = TRUE,
#   print_plots = TRUE
# )
# lthal_lpme_laplacian_est <- calc_lpme_est_params(lthal_lpme_result_laplacian, lthal_mat)
# lthal_lpme_laplacian_vals <- lthal_lpme_laplacian_est$results
# lthal_lpme_laplacian_params <- lthal_lpme_laplacian_est$params

lthal_pme_result <- list()
lthal_pme_vals <- list()
for (t in 1:length(time_vals)) {
  temp_data <- lthal_mat[lthal_mat[, 1] == time_vals[t], -1]
  lthal_pme_result[[t]] <- pme(temp_data, d = 2, lambda = exp(-12:5), verbose = TRUE)
  lthal_pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(lthal_pme_result[[t]], temp_data))
}
lthal_pme_vals <- reduce(lthal_pme_vals, rbind)

write.csv(lthal_mat, "results/adni_lthal_mat.csv")
write.csv(lthal_lpme_isomap_vals, "results/adni_lthal_lpme_isomap_vals.csv")
write.csv(lthal_pme_vals, "results/adni_lthal_pme_vals.csv")

p <- plot_ly(
  x = lthal_mat[, 2],
  y = lthal_mat[, 3],
  z = lthal_mat[, 4],
  frame = lthal_mat[, 1],
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3),
  name = "Data"
) %>%
  add_markers(
    x = lthal_lpme_isomap_vals[, 2],
    y = lthal_lpme_isomap_vals[, 3],
    z = lthal_lpme_isomap_vals[, 4],
    frame = lthal_lpme_isomap_vals[, 1],
    name = "LPME"
  ) %>%
  add_markers(
    x = lthal_pme_vals[, 2],
    y = lthal_pme_vals[, 3],
    z = lthal_pme_vals[, 4],
    frame = lthal_pme_vals[, 1],
    name = "PME"
  )

lthal_lpme_plots <- list()
lthal_pme_plots <- list()
lthal_data_plots <- list()

colors <- brewer.pal(3, "Set1")

png(
  "results/adni_plots/adni_lthal_data_plot.png",
  res = 1000,
  height = 15000,
  width = 3500
)
par(oma = c(4, 1, 1, 1), mfrow = c(6, 1), mar = c(2, 2, 1, 1))

for (t in 1:length(time_vals)) {
  temp_data <- lthal_mat[lthal_mat[, 1] == time_vals[t], ]
  scatter3D(
    x = temp_data[, 2],
    y = temp_data[, 3],
    z = temp_data[, 4],
    pch = 20,
    col = alpha(colors[1], 0.5),
    theta = 220,
    ticktype = "detailed",
    main = paste0("Time = ", round((time_vals * max(lthal$duration))[t], 2)),
    xlim = c(-0.15, 0.15),
    ylim = c(-0.15, 0.15),
    zlim = c(-0.15, 0.15)
  )
}
dev.off()

png(
  "results/adni_plots/adni_lthal_lpme_isomap_plot.png",
  res = 1000,
  height = 15000,
  width = 3500
)
par(oma = c(4, 1, 1, 1), mfrow = c(6, 1), mar = c(2, 2, 1, 1))

for (t in 1:length(time_vals)) {
  temp_lpme <- lthal_lpme_isomap_vals[lthal_lpme_isomap_vals[, 1] == time_vals[t], ]
  scatter3D(
    x = temp_lpme[, 2],
    y = temp_lpme[, 3],
    z = temp_lpme[, 4],
    pch = 20,
    col = alpha(colors[2], 0.5),
    theta = 220,
    ticktype = "detailed",
    main = paste0("Time = ", round((time_vals * max(lthal$duration))[t], 2)),
    xlim = c(-0.15, 0.15),
    ylim = c(-0.15, 0.15),
    zlim = c(-0.15, 0.15)
  )
}
dev.off()

png(
  "results/adni_plots/adni_lthal_lpme_dm_plot.png",
  res = 1000,
  height = 15000,
  width = 3500
)
par(oma = c(4, 1, 1, 1), mfrow = c(6, 1), mar = c(2, 2, 1, 1))

for (t in 1:length(time_vals)) {
  temp_lpme <- lthal_lpme_dm_vals[lthal_lpme_dm_vals[, 1] == time_vals[t], ]
  scatter3D(
    x = temp_lpme[, 2],
    y = temp_lpme[, 3],
    z = temp_lpme[, 4],
    pch = 20,
    col = alpha(colors[2], 0.5),
    theta = 220,
    ticktype = "detailed",
    main = paste0("Time = ", round((time_vals * max(lthal$duration))[t], 2)),
    xlim = c(-0.15, 0.15),
    ylim = c(-0.15, 0.15),
    zlim = c(-0.15, 0.15)
  )
}
dev.off()

png(
  "results/adni_plots/adni_lthal_lpme_laplacian_plot.png",
  res = 1000,
  height = 15000,
  width = 3500
)
par(oma = c(4, 1, 1, 1), mfrow = c(6, 1), mar = c(2, 2, 1, 1))

for (t in 1:length(time_vals)) {
  temp_lpme <- lthal_lpme_laplacian_vals[lthal_lpme_laplacian_vals[, 1] == time_vals[t], ]
  scatter3D(
    x = temp_lpme[, 2],
    y = temp_lpme[, 3],
    z = temp_lpme[, 4],
    pch = 20,
    col = alpha(colors[2], 0.5),
    theta = 220,
    ticktype = "detailed",
    main = paste0("Time = ", round((time_vals * max(lthal$duration))[t], 2)),
    xlim = c(-0.15, 0.15),
    ylim = c(-0.15, 0.15),
    zlim = c(-0.15, 0.15)
  )
}
dev.off()

png(
  "results/adni_plots/adni_lthal_pme_plot.png",
  res = 1000,
  height = 15000,
  width = 3500
)
par(oma = c(4, 1, 1, 1), mfrow = c(6, 1), mar = c(2, 2, 1, 1))

for (t in 1:length(time_vals)) {
  temp_pme <- lthal_pme_vals[lthal_pme_vals[, 1] == time_vals[t], ]
  scatter3D(
    x = temp_pme[, 2],
    y = temp_pme[, 3],
    z = temp_pme[, 4],
    pch = 20,
    col = alpha(colors[3], 0.5),
    theta = 220,
    ticktype = "detailed",
    main = paste0("Time = ", round((time_vals * max(lthal$duration))[t], 2)),
    xlim = c(-0.15, 0.15),
    ylim = c(-0.15, 0.15),
    zlim = c(-0.15, 0.15)
  )
}
dev.off()

lthal_cross_section <- create_cross_section_matrix(
  list(lthal_mat, lthal_lpme_isomap_vals, lthal_pme_vals),
  2,
  c(0.005, 0.005, 0.005),
  nrow = 2,
  lhipp = lthal
)
ggsave(
  "results/adni_plots/adni_lthal_cross_section.png",
  plot = lthal_cross_section,
  height = 5000,
  width = 7000,
  units = "px",
  dpi = 1250
)

# #### VOLUME ESTIMATION
#
# # Voxel dimension: 1.2mm x 0.9375mm x 0.9375mm
# # Voxel volume: 1.054688 mm^3
#
# # Begin by rescaling data and estimated surfaces to full size
#
# lthal_data_rescaled <- list()
# lthal_lpme_rescaled <- list()
# lthal_pme_rescaled <- list()
#
# for (time_idx in 1:length(time_vals)) {
#   temp_data <- lthal_mat[lthal_mat[, 1] == time_vals[time_idx], ]
#   temp_lpme <- lthal_lpme_vals[lthal_lpme_vals[, 1] == time_vals[time_idx], ]
#   temp_pme <- lthal_pme_vals[lthal_pme_vals[, 1] == time_vals[time_idx], ]
#
#   temp_max_x <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$max_x)
#   temp_max_y <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$max_y)
#   temp_max_z <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$max_z)
#   # temp_sd_x <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$sd_x)
#   # temp_sd_y <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$sd_y)
#   # temp_sd_z <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$sd_z)
#
#   temp_data_rescaled <- temp_data[, -(5:7)]
#   temp_lpme_rescaled <- temp_lpme[, -(5:7)]
#   temp_pme_rescaled <- temp_pme[, -(5:7)]
#
#   temp_data_rescaled[, 2] <- round((temp_data_rescaled[, 2] * temp_max_x) / 1.2)
#   temp_data_rescaled[, 3] <- round((temp_data_rescaled[, 3] * temp_max_y) / 0.9375)
#   temp_data_rescaled[, 4] <- round((temp_data_rescaled[, 4] * temp_max_z) / 0.9375)
#
#   temp_lpme_rescaled[, 2] <- round((temp_lpme_rescaled[, 2] * temp_max_x) / 1.2)
#   temp_lpme_rescaled[, 3] <- round((temp_lpme_rescaled[, 3] * temp_max_y) / 0.9375)
#   temp_lpme_rescaled[, 4] <- round((temp_lpme_rescaled[, 4] * temp_max_z) / 0.9375)
#
#   temp_pme_rescaled[, 2] <- round((temp_pme_rescaled[, 2] * temp_max_x) / 1.2)
#   temp_pme_rescaled[, 3] <- round((temp_pme_rescaled[, 3] * temp_max_y) / 0.9375)
#   temp_pme_rescaled[, 4] <- round((temp_pme_rescaled[, 4] * temp_max_z) / 0.9375)
#
#   # temp_data_rescaled[, 2] <- round((temp_data_rescaled[, 2] * temp_sd_x) / 1.2)
#   # temp_data_rescaled[, 3] <- round((temp_data_rescaled[, 3] * temp_sd_y) / 0.9375)
#   # temp_data_rescaled[, 4] <- round((temp_data_rescaled[, 4] * temp_sd_z) / 0.9375)
#   #
#   # temp_lpme_rescaled[, 2] <- round((temp_lpme_rescaled[, 2] * temp_sd_x) / 1.2)
#   # temp_lpme_rescaled[, 3] <- round((temp_lpme_rescaled[, 3] * temp_sd_y) / 0.9375)
#   # temp_lpme_rescaled[, 4] <- round((temp_lpme_rescaled[, 4] * temp_sd_z) / 0.9375)
#   #
#   # temp_pme_rescaled[, 2] <- round((temp_pme_rescaled[, 2] * temp_sd_x) / 1.2)
#   # temp_pme_rescaled[, 3] <- round((temp_pme_rescaled[, 3] * temp_sd_y) / 0.9375)
#   # temp_pme_rescaled[, 4] <- round((temp_pme_rescaled[, 4] * temp_sd_z) / 0.9375)
#
#   lthal_data_rescaled[[time_idx]] <- temp_data_rescaled
#   lthal_lpme_rescaled[[time_idx]] <- temp_lpme_rescaled
#   lthal_pme_rescaled[[time_idx]] <- temp_pme_rescaled
# }
#
# lthal_data_rescaled_full <- reduce(lthal_data_rescaled, rbind)
# lthal_lpme_rescaled_full <- reduce(lthal_lpme_rescaled, rbind)
# lthal_pme_rescaled_full <- reduce(lthal_pme_rescaled, rbind)
#
# voxel_vol <- 1.2 * 0.9375 * 0.9375
#
# lthal_est_data_volume <- map(
#   lthal_data_rescaled,
#   ~ estimate_volume(.x[, -1], voxel_vol)
# )
#
# lthal_est_data_volume_points <- map(
#   lthal_est_data_volume,
#   ~ .x[[2]]
# )
#
# lthal_est_data_volume <- map(
#   lthal_est_data_volume,
#   ~ .x[[1]]
# ) %>%
#   reduce(c)
#
# lthal_est_lpme_volume <- map(
#   lthal_lpme_rescaled,
#   ~ estimate_volume(.x[, -1], voxel_vol)
# )
#
# lthal_est_lpme_volume_points <- map(
#   lthal_est_lpme_volume,
#   ~ .x[[2]]
# )
#
# lthal_est_lpme_volume <- map(
#   lthal_est_lpme_volume,
#   ~ .x[[1]]
# ) %>%
#   reduce(c)
#
# lthal_est_pme_volume <- map(
#   lthal_pme_rescaled,
#   ~ estimate_volume(.x[, -1], voxel_vol)
# )
#
# lthal_est_pme_volume_points <- map(
#   lthal_est_pme_volume,
#   ~ .x[[2]]
# )
#
# lthal_est_pme_volume <- map(
#   lthal_est_pme_volume,
#   ~ .x[[1]]
# ) %>%
#   reduce(c)
#
# lthal_est_data_volume_points_full <- map(
#   1:length(lthal_est_data_volume_points),
#   ~ cbind(time_vals[.x], lthal_est_data_volume_points[[.x]])
# ) %>%
#   reduce(rbind)
#
#
# lthal_est_lpme_volume_points_full <- map(
#   1:length(lthal_est_lpme_volume_points),
#   ~ cbind(time_vals[.x], lthal_est_lpme_volume_points[[.x]])
# ) %>%
#   reduce(rbind)
#
# lthal_est_pme_volume_points_full <- map(
#   1:length(lthal_est_pme_volume_points),
#   ~ cbind(time_vals[.x], lthal_est_pme_volume_points[[.x]])
# ) %>%
#   reduce(rbind)
#

lthal_rescaled <- lthal %>%
  mutate(
    x_rescaled = round((x * max_x)),
    y_rescaled = round((y * max_y)),
    z_rescaled = round((z * max_z))
  )
candidate_x <- seq(min(lthal_rescaled$x_rescaled), max(lthal_rescaled$x_rescaled), 1)
candidate_y <- seq(min(lthal_rescaled$y_rescaled), max(lthal_rescaled$y_rescaled), 1)
candidate_z <- seq(min(lthal_rescaled$z_rescaled), max(lthal_rescaled$z_rescaled), 1)

candidate_voxels <- expand_grid(candidate_x, candidate_y, candidate_z)

lthal_pme_volumes <- vector(mode = "numeric", length = length(time_vals))
lthal_lpme_volumes <- vector(mode = "numeric", length = length(time_vals))
lthal_pme_interior_plots <- list()
lthal_lpme_interior_plots <- list()

for (time_idx in 1:length(time_vals)) {
  temp_max_x <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$max_x)
  temp_max_y <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$max_y)
  temp_max_z <- unique(lthal[lthal$time_from_bl == time_vals[time_idx], ]$max_z)

  temp_candidates <- candidate_voxels %>%
    mutate(
      x_scaled = candidate_x / temp_max_x,
      y_scaled = candidate_y / temp_max_y,
      z_scaled = candidate_z / temp_max_z
    )

  temp_pme <- lthal_pme_result[[time_idx]]

  temp_pme_embedding <- temp_pme$embedding_map
  temp_pme_coefs <- temp_pme$coefs[[which.min(temp_pme$MSD)]]
  temp_pme_params <- temp_pme$parameterization[[which.min(temp_pme$MSD)]]

  temp_lpme_coefs <- lthal_lpme_result_isomap$sol_coef_functions[[which.min(lthal_lpme_result_isomap$msd)]](time_vals[time_idx])
  temp_lpme_params <- lthal_lpme_result_isomap$parameterization_list[[which.min(lthal_lpme_result_isomap$msd)]]
  temp_lpme_params <- temp_lpme_params[temp_lpme_params[, 1] == time_vals[time_idx], -1]
  lpme_n_knots <- nrow(temp_lpme_params)
  d <- ncol(temp_lpme_params)
  coef_mat <- matrix(
    temp_lpme_coefs,
    lpme_n_knots + d + 1,
    byrow = TRUE
  )

  temp_lpme_embedding <- function(r) {
    t(coef_mat[1:lpme_n_knots, ]) %*% pme::etaFunc(r, temp_lpme_params, 4 - d) +
      t(coef_mat[(lpme_n_knots + 1):(lpme_n_knots + d + 1), ]) %*% matrix(c(1, r), ncol = 1)
  }

  lthal_interior_voxel_pme <- vector(length = nrow(temp_candidates))
  lthal_interior_voxel_lpme <- vector(length = nrow(temp_candidates))
  for (row_idx in 1:nrow(temp_candidates)) {
    lthal_interior_voxel_pme[row_idx] <- interior_identification(
      temp_pme_embedding,
      temp_pme_coefs,
      temp_pme_params,
      unlist(temp_candidates[row_idx, 4:6]),
      c(0, 0, 0)
    )
    lthal_interior_voxel_lpme[row_idx] <- interior_identification(
      temp_lpme_embedding,
      coef_mat,
      temp_lpme_params,
      unlist(temp_candidates[row_idx, 4:6]),
      c(0, 0, 0)
    )
  }

  lthal_interior_points_pme <- temp_candidates[lthal_interior_voxel_pme, ]
  lthal_interior_points_lpme <- temp_candidates[lthal_interior_voxel_lpme, ]
  p <- plot_ly(
    lthal_interior_points_pme,
    x = ~candidate_x,
    y = ~candidate_y,
    z = ~candidate_z,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  )

  p_lpme <- plot_ly(
    lthal_interior_points_lpme,
    x = ~candidate_x,
    y = ~candidate_y,
    z = ~candidate_z,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3)
  )

  lthal_pme_interior_plots[[time_idx]] <- p
  lthal_lpme_interior_plots[[time_idx]] <- p_lpme
  lthal_pme_volumes[time_idx] <- 1.054688 * nrow(lthal_interior_points_pme)
  lthal_lpme_volumes[time_idx] <- 1.054688 * nrow(lthal_interior_points_lpme)
}

adni_hipp_info <- read_csv("data/adni_fsl_hipp_info.csv")
adni_thal_info <- read_csv("data/adni_fsl_thal_red_info.csv")

lhipp_est_data_volume <- adni_hipp_info %>%
  filter(patno == "002_S_0685") %>%
  select(lhipp_vol) %>%
  rename(data = lhipp_vol)

lthal_est_data_volume <- adni_thal_info %>%
  filter(patno == "002_S_0685") %>%
  select(lthal_vol) %>%
  rename(data = lthal_vol)

lhipp_volume_df <- data.frame(
  time = time_vals,
  data = lhipp_est_data_volume,
  lpme = lhipp_lpme_volumes,
  pme = lhipp_pme_volumes
)

write_csv(lhipp_volume_df, "results/lhipp_volume.csv")

lthal_volume_df <- data.frame(
  time = time_vals,
  data = lthal_est_data_volume,
  lpme = lthal_lpme_volumes,
  pme = lthal_pme_volumes
)

write_csv(lthal_volume_df, "results/lthal_volume.csv")

lhipp_volume <- lhipp_volume_df %>%
  pivot_longer(cols = c(data, lpme, pme), names_to = "Source", values_to = "Volume") %>%
  rename(Time = time) %>%
  mutate(Source = str_to_title(Source))

lthal_volume <- lthal_volume_df %>%
  pivot_longer(cols = c(data, lpme, pme), names_to = "Source", values_to = "Volume") %>%
  rename(Time = time) %>%
  mutate(Source = str_to_title(Source))

ggplot(lhipp_volume) +
  geom_point(aes(x = Time, y = Volume, color = Source)) +
  geom_line(aes(x = Time, y = Volume, color = Source)) +
  ylab("Estimated Volume") +
  scale_color_manual(values = colors)
ggsave("results/adni_plots/lhipp_volume.png")

ggplot(lthal_volume) +
  geom_point(aes(x = Time, y = Volume, color = Source)) +
  geom_line(aes(x = Time, y = Volume, color = Source)) +
  ylab("Estimated Volume") +
  scale_color_manual(values = colors)
ggsave("results/adni_plots/lthal_volume.png")
