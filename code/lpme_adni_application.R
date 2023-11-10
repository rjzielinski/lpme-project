library(lubridate)
library(plotly)
library(plot3D)
library(pme)
library(pracma)
library(RColorBrewer)
library(Rfast)
library(tidyverse)

source("code/functions/calc_pme_est.R")
source("code/functions/calc_lpme_est.R")

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

lhipp_surface <- read_csv("data/adni_fsl_lhipp_surface_red.csv")

# process scan dates
lhipp_surface <- lhipp_surface %>%
  mutate(
    scan_date = gsub("\\.0", "", scan_date),
    scan_date = ymd_hms(scan_date),
    scan_date = decimal_date(scan_date)
  )

lhipp_surface_centers <- lhipp_surface %>%
  group_by(patno, scan_date) %>%
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
    time_from_bl = (scan_date - time_bl) / duration
  )

lhipp_surface_spherical <- lhipp_surface %>%
  dplyr::select(x, y, z) %>%
  as.matrix() %>%
  cart2sph() %>%
  as_tibble() %>%
  mutate(
    theta = theta / max(abs(theta)),
    phi = phi / max(abs(phi)),
    r = r / max(abs(r))
  )

lhipp_surface_full <- bind_cols(lhipp_surface, lhipp_surface_spherical)

lhipp <- lhipp_surface_full %>%
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

time_vals <- unique(lhipp_mat[, 1])

set.seed(1036)
lpme_result <- lpme(lhipp_mat, d = 2, verbose = TRUE, print_plots = TRUE)
lpme_vals <- calc_lpme_est(lpme_result, lhipp_mat)

pme_result <- list()
pme_vals <- list()
for (t in 1:length(time_vals)) {
  temp_data <- lhipp_mat[lhipp_mat[, 1] == time_vals[t], -1]
  pme_result[[t]] <- pme(temp_data, d = 2, verbose = "MSD")
  pme_vals[[t]] <- cbind(time_vals[t], calc_pme_est(pme_result[[t]], temp_data))
}
pme_vals <- reduce(pme_vals, rbind)

write.csv(lhipp_mat, "results/adni_lhipp_mat.csv")
write.csv(lpme_vals, "results/adni_lpme_vals.csv")
write.csv(pme_vals, "results/adni_pme_vals.csv")

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
    x = lpme_vals[, 2],
    y = lpme_vals[, 3],
    z = lpme_vals[, 4],
    frame = lpme_vals[, 1],
    name = "LPME"
  ) %>%
  add_markers(
    x = pme_vals[, 2],
    y = pme_vals[, 3],
    z = pme_vals[, 4],
    frame = pme_vals[, 1],
    name = "PME"
  )

lpme_plots <- list()
pme_plots <- list()
data_plots <- list()

colors <- brewer.pal(3, "Set1")

png(
  "results/adni_plots/adni_data_plot.png",
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
  "results/adni_plots/adni_lpme_plot.png",
  res = 1000,
  height = 15000,
  width = 3500
)
par(oma = c(4, 1, 1, 1), mfrow = c(6, 1), mar = c(2, 2, 1, 1))

for (t in 1:length(time_vals)) {
  temp_lpme <- lpme_vals[lpme_vals[, 1] == time_vals[t], ]
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
  "results/adni_plots/adni_pme_plot.png",
  res = 1000,
  height = 15000,
  width = 3500
)
par(oma = c(4, 1, 1, 1), mfrow = c(6, 1), mar = c(2, 2, 1, 1))

for (t in 1:length(time_vals)) {
  temp_pme <- pme_vals[pme_vals[, 1] == time_vals[t], ]
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
  list(lhipp_mat, lpme_vals, pme_vals),
  2,
  c(0.005, 0.005, 0.005),
  nrow = 2,
  lhipp = lhipp
)
ggsave(
  "results/adni_plots/adni_cross_section.png",
  plot = lhipp_cross_section,
  height = 5000,
  width = 7000,
  units = "px",
  dpi = 1250
)
