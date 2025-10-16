library(brolgar)
library(dplyr)
library(here)
library(plotly)
library(pme)
library(readr)
library(reticulate)
library(RColorBrewer)

use_condaenv("lpme")
np <- import("numpy")
o3d <- import("open3d")
pv <- import("pyvista")
# pymeshfix <- import("pymeshfix")

# source(here("code/functions/estimate_mesh_volume.R"))
source(here("code/functions/estimate_mesh_volume_poisson.R"))

result_dir <- here("output/adni")

isolate_patno <- function(file_vec, file_path, structure) {
  patnos <- gsub(file_path, "", file_vec)
  patnos <- gsub(paste0(structure, "_results.RDS"), "", patnos)
  patnos <- gsub("/", "", patnos)
}

result_files <- list.files(result_dir, full.names = TRUE, recursive = TRUE)

lhipp_files <- result_files[grepl("lhipp", result_files)]
rhipp_files <- result_files[grepl("rhipp", result_files)]
lthal_files <- result_files[grepl("lthal", result_files)]
rthal_files <- result_files[grepl("rthal", result_files)]

lhipp_patnos <- isolate_patno(lhipp_files, result_dir, "lhipp")
rhipp_patnos <- isolate_patno(rhipp_files, result_dir, "rhipp")
lthal_patnos <- isolate_patno(lthal_files, result_dir, "lthal")
rthal_patnos <- isolate_patno(rthal_files, result_dir, "rthal")

# common_patnos <- intersect(lhipp_patnos, lthal_patnos)
common_patnos <- lhipp_patnos

lhipp_idx <- which(grepl(common_patnos[4], lhipp_files))
lthal_idx <- which(grepl(common_patnos[4], lthal_files))

lhipp_ex <- readRDS(lhipp_files[lhipp_idx])
lthal_ex <- readRDS(lthal_files[lthal_idx])

structure_plot_grid <- function(adni_est, silhouette_decimation = 0.75) {
  time_points <- unique(adni_est$data$time_from_bl)

  palette <- colorRampPalette(c(brewer.pal(8, "Spectral")))
  palette_colors <- palette(101)
  time_colors <- palette_colors[ceiling((time_points + 1e-10) * 100)]

  adni_plot <- pv$Plotter(shape = tuple(4L, as.integer(length(time_points))))
  adni_sil_plot <- pv$Plotter(shape = tuple(2L, 2L))

  for (time_idx in seq_along(time_points)) {
    temp_data <- adni_est$data |>
      filter(time_from_bl == time_points[time_idx]) |>
      select(x, y, z) |>
      as.matrix()

    temp_lpme_part_recon <- adni_est$lpme_part$reconstructions[
      adni_est$lpme_part$reconstructions[, 1] == time_points[time_idx],
      -1
    ]
    temp_pme_part_recon <- adni_est$pme_part$reconstructions[
      adni_est$pme_part$reconstructions[, 1] == time_points[time_idx],
      -1
    ]
    temp_pc_part_recon <- adni_est$pc_part$reconstructions[
      adni_est$pc_part$reconstructions[, 1] == time_points[time_idx],
      -1
    ]

    data_mesh <- estimate_mesh_volume_poisson(temp_data)
    lpme_mesh <- estimate_mesh_volume_poisson(temp_lpme_part_recon)
    pme_mesh <- estimate_mesh_volume_poisson(temp_pme_part_recon)
    pc_mesh <- estimate_mesh_volume_poisson(temp_pc_part_recon)

    adni_plot$subplot(0L, as.integer(time_idx - 1))
    adni_plot$add_mesh(
      # data_mesh$mesh$decimate(0.9, preserve_topology = TRUE),
      data_mesh$mesh,
      color = "gray",
      show_edges = TRUE,
      opacity = 0.5
    )

    adni_plot$add_text(
      paste("Time =", round(time_points[time_idx], 2)),
      position = "upper_left",
      font_size = 12
    )

    if (time_idx == 1) {
      adni_plot$add_text(
        "Data",
        position = "lower_left",
        font_size = 12
      )
    }

    adni_sil_plot$subplot(0L, 0L)
    adni_sil_plot$add_silhouette(
      data_mesh$mesh,
      color = time_colors[time_idx],
      decimate = silhouette_decimation
    )

    if (time_idx == 1) {
      adni_sil_plot$add_text(
        "Data",
        position = "lower_left",
        font_size = 12
      )
    }

    adni_plot$subplot(1L, as.integer(time_idx - 1))
    adni_plot$add_mesh(
      lpme_mesh$mesh,
      color = "blue",
      show_edges = TRUE,
      opacity = 0.75
    )

    if (time_idx == 1) {
      adni_plot$add_text(
        "LPME",
        position = "lower_left",
        font_size = 12
      )
    }

    adni_sil_plot$subplot(0L, 1L)
    adni_sil_plot$add_silhouette(
      lpme_mesh$mesh,
      color = time_colors[time_idx],
      decimate = silhouette_decimation
    )

    if (time_idx == 1) {
      adni_sil_plot$add_text(
        "LPME",
        position = "lower_left",
        font_size = 12
      )
    }

    adni_plot$subplot(2L, as.integer(time_idx - 1))
    adni_plot$add_mesh(
      pme_mesh$mesh,
      color = "green",
      show_edges = TRUE,
      opacity = 0.75
    )

    if (time_idx == 1) {
      adni_plot$add_text(
        "PME",
        position = "lower_left",
        font_size = 12
      )
    }

    adni_sil_plot$subplot(1L, 0L)
    adni_sil_plot$add_silhouette(
      pme_mesh$mesh,
      color = time_colors[time_idx],
      decimate = silhouette_decimation
    )

    if (time_idx == 1) {
      adni_sil_plot$add_text(
        "PME",
        position = "lower_left",
        font_size = 12
      )
    }

    adni_plot$subplot(3L, as.integer(time_idx - 1))
    adni_plot$add_mesh(
      pc_mesh$mesh,
      color = "purple",
      show_edges = TRUE,
      opacity = 0.75
    )

    if (time_idx == 1) {
      adni_plot$add_text(
        "PS",
        position = "lower_left",
        font_size = 12
      )
    }

    adni_sil_plot$subplot(1L, 1L)
    adni_sil_plot$add_silhouette(
      pc_mesh$mesh,
      color = time_colors[time_idx],
      decimate = silhouette_decimation
    )

    if (time_idx == 1) {
      adni_sil_plot$add_text(
        "PS",
        position = "lower_left",
        font_size = 12
      )
    }
  }
  adni_plot$link_views()
  adni_plot$camera_position <- "yz"
  adni_plot$camera$azimuth <- -45

  adni_sil_plot$link_views()
  adni_sil_plot$camera_position <- "yz"
  adni_sil_plot$camera$azimuth <- -45

  list(mesh_plot = adni_plot, sil_plot = adni_sil_plot)
}

lhipp_plots <- structure_plot_grid(lhipp_ex)

lhipp_plots$mesh_plot$show()
lhipp_plots$mesh_plot$screenshot("test_lhipp_mesh_plot.png")

lhipp_plots$sil_plot$show()
lhipp_plots$sil_plot$screenshot("test_lhipp_sil_plot.png")

lthal_plots <- structure_plot_grid(lthal_ex)

lthal_plots$mesh_plot$show()
lthal_plots$mesh_plot$screenshot("test_lthal_mesh_plot.png")

colors <- brewer.pal(3, "Set1")


# RETRIEVE VOLUME ESTIMATES FROM MESHES -------------------------------------

adni <- read_csv("data/adni_info_full.csv")

lhipp_surface <- read_csv(here("data/lhipp_surface_fsl_processed.csv"))
rhipp_surface <- read_csv(here("data/rhipp_surface_fsl_processed.csv"))
lthal_surface <- read_csv(here("data/lthal_surface_fsl_processed.csv"))
rthal_surface <- read_csv(here("data/rthal_surface_fsl_processed.csv"))

# adni_status <- full_join(
#   adni_status,
#   roster,
#   by = c("Phase", "ID", "RID", "SITEID")
# )

hipp_info <- data.frame(
  patno = character(),
  date = numeric(),
  lhipp_data_vol = numeric(),
  lhipp_lpme_vol = numeric(),
  lhipp_pme_vol = numeric(),
  lhipp_pc_vol = numeric(),
  rhipp_data_vol = numeric(),
  rhipp_lpme_vol = numeric(),
  rhipp_pme_vol = numeric(),
  rhipp_pc_vol = numeric()
)

thal_info <- data.frame(
  patno = character(),
  date = numeric(),
  lthal_data_vol = numeric(),
  lthal_lpme_vol = numeric(),
  lthal_pme_vol = numeric(),
  lthal_pc_vol = numeric(),
  rthal_data_vol = numeric(),
  rthal_lpme_vol = numeric(),
  rthal_pme_vol = numeric(),
  rthal_pc_vol = numeric()
)


for (patno in common_patnos) {
  print(
    paste0(
      "PATNO ",
      which(common_patnos == patno),
      " of ",
      length(common_patnos),
      ": ",
      patno
    )
  )
  patno_dir <- paste0(result_dir, "/", patno, "/")
  lhipp_info <- lhipp_surface |>
    filter(subid == patno)
  rhipp_info <- rhipp_surface |>
    filter(subid == patno)
  lthal_info <- lthal_surface |>
    filter(subid == patno)
  rthal_info <- rthal_surface |>
    filter(subid == patno)

  patno_dates <- unique(c(
    lhipp_info$date,
    rhipp_info$date,
    lthal_info$date,
    rthal_info$date
  ))

  lhipp_result <- readRDS(paste0(patno_dir, "lhipp_results.RDS"))
  rhipp_result <- readRDS(paste0(patno_dir, "rhipp_results.RDS"))
  lthal_result <- readRDS(paste0(patno_dir, "lthal_results.RDS"))
  rthal_result <- readRDS(paste0(patno_dir, "rthal_results.RDS"))

  time_points <- unique(lhipp_result$data$time_from_bl)
  lhipp_data_vol <- vector(mode = "numeric", length = length(time_points))
  rhipp_data_vol <- vector(mode = "numeric", length = length(time_points))
  lthal_data_vol <- vector(mode = "numeric", length = length(time_points))
  rthal_data_vol <- vector(mode = "numeric", length = length(time_points))

  for (time_idx in seq_along(time_points)) {
    lhipp_scaling <- lhipp_info |>
      filter(time_from_bl == time_points[time_idx]) |>
      summarize(max_x = mean(max_x), max_y = mean(max_y), max_z = mean(max_z))
    rhipp_scaling <- rhipp_info |>
      filter(time_from_bl == time_points[time_idx]) |>
      summarize(max_x = mean(max_x), max_y = mean(max_y), max_z = mean(max_z))
    lthal_scaling <- lthal_info |>
      filter(time_from_bl == time_points[time_idx]) |>
      summarize(max_x = mean(max_x), max_y = mean(max_y), max_z = mean(max_z))
    rthal_scaling <- rthal_info |>
      filter(time_from_bl == time_points[time_idx]) |>
      summarize(max_x = mean(max_x), max_y = mean(max_y), max_z = mean(max_z))

    temp_lhipp_data <- lhipp_result$data |>
      filter(time_from_bl == time_points[time_idx]) |>
      select(x, y, z) |>
      as.matrix()
    temp_rhipp_data <- rhipp_result$data |>
      filter(time_from_bl == time_points[time_idx]) |>
      select(x, y, z) |>
      as.matrix()

    temp_lthal_data <- lthal_result$data |>
      filter(time_from_bl == time_points[time_idx]) |>
      select(x, y, z) |>
      as.matrix()
    temp_rthal_data <- rthal_result$data |>
      filter(time_from_bl == time_points[time_idx]) |>
      select(x, y, z) |>
      as.matrix()

    lhipp_data_vol[time_idx] <- estimate_mesh_volume_poisson(
      temp_lhipp_data
    )$volume *
      (lhipp_scaling$max_x *
        lhipp_scaling$max_y *
        lhipp_scaling$max_z)
    rhipp_data_vol[time_idx] <- estimate_mesh_volume_poisson(
      temp_rhipp_data
    )$volume *
      (rhipp_scaling$max_x *
        rhipp_scaling$max_y *
        rhipp_scaling$max_z)

    lthal_data_vol[time_idx] <- estimate_mesh_volume_poisson(
      temp_lthal_data
    )$volume *
      (lthal_scaling$max_x *
        lthal_scaling$max_y *
        lthal_scaling$max_z)
    rthal_data_vol[time_idx] <- estimate_mesh_volume_poisson(
      temp_rthal_data
    )$volume *
      (rthal_scaling$max_x *
        rthal_scaling$max_y *
        rthal_scaling$max_z)
  }

  lhipp_lpme_vols <- lhipp_result$lpme_part$volumes_mesh
  lhipp_pme_vols <- lhipp_result$pme_part$volumes_mesh
  lhipp_pc_vols <- lhipp_result$pc_part$volumes

  rhipp_lpme_vols <- rhipp_result$lpme_part$volumes_mesh
  rhipp_pme_vols <- rhipp_result$pme_part$volumes_mesh
  rhipp_pc_vols <- rhipp_result$pc_part$volumes

  lthal_lpme_vols <- lthal_result$lpme_part$volumes_mesh
  lthal_pme_vols <- lthal_result$pme_part$volumes_mesh
  lthal_pc_vols <- lthal_result$pc_part$volumes

  rthal_lpme_vols <- rthal_result$lpme_part$volumes_mesh
  rthal_pme_vols <- rthal_result$pme_part$volumes_mesh
  rthal_pc_vols <- rthal_result$pc_part$volumes

  patno_hipp_info <- data.frame(
    patno = rep(patno, length(patno_dates)),
    date = patno_dates,
    lhipp_data_vol = lhipp_data_vol,
    lhipp_lpme_vol = lhipp_lpme_vols,
    lhipp_pme_vol = lhipp_pme_vols,
    lhipp_pc_vol = lhipp_pc_vols,
    rhipp_data_vol = rhipp_data_vol,
    rhipp_lpme_vol = rhipp_lpme_vols,
    rhipp_pme_vol = rhipp_pme_vols,
    rhipp_pc_vol = rhipp_pc_vols
  )
  patno_thal_info <- data.frame(
    patno = rep(patno, length(patno_dates)),
    date = patno_dates,
    lthal_data_vol = lthal_data_vol,
    lthal_lpme_vol = lthal_lpme_vols,
    lthal_pme_vol = lthal_pme_vols,
    lthal_pc_vol = lthal_pc_vols,
    rthal_data_vol = rthal_data_vol,
    rthal_lpme_vol = rthal_lpme_vols,
    rthal_pme_vol = rthal_pme_vols,
    rthal_pc_vol = rthal_pc_vols
  )

  hipp_info <- rbind(hipp_info, patno_hipp_info)
  thal_info <- rbind(thal_info, patno_thal_info)
}

write_csv(hipp_info, here("output/test_adni_hipp_volumes.csv"))
write_csv(thal_info, here("output/test_adni_thal_volumes.csv"))


hipp_bl <- hipp_info %>%
  group_by(patno) %>%
  summarize(
    date_bl = min(date),
    max_date = max(date)
  )
thal_bl <- thal_info %>%
  group_by(patno) %>%
  summarize(
    date_bl = min(date),
    max_date = max(date)
  )

eligible_patnos <- hipp_bl %>%
  filter(max_date - date_bl > 2) %>%
  select(patno) %>%
  unlist() %>%
  unique()

hipp_info <- full_join(hipp_info, hipp_bl, by = "patno") %>%
  mutate(time_from_bl = date - date_bl) %>%
  filter(patno %in% eligible_patnos)

thal_info <- full_join(thal_info, thal_bl, by = "patno") %>%
  mutate(time_from_bl = date - date_bl) %>%
  filter(patno %in% eligible_patnos)

hipp_summary <- hipp_info %>%
  group_by(patno) %>%
  summarize(
    lhipp_data_vol2_sd = sd(lhipp_data_vol2),
    lhipp_vol_lpme1_sd = sd(lhipp_vol_lpme1),
    lhipp_vol_lpme2_sd = sd(lhipp_vol_lpme2),
    lhipp_vol_pme1_sd = sd(lhipp_vol_pme1),
    lhipp_vol_pme2_sd = sd(lhipp_vol_pme2),
    rhipp_data_vol2_sd = sd(rhipp_data_vol2),
    rhipp_vol_lpme1_sd = sd(rhipp_vol_lpme1),
    rhipp_vol_lpme2_sd = sd(rhipp_vol_lpme2),
    rhipp_vol_pme1_sd = sd(rhipp_vol_pme1),
    rhipp_vol_pme2_sd = sd(rhipp_vol_pme2)
  )

hipp_long <- hipp_info %>%
  pivot_longer(
    lhipp_data_vol2:lhipp_vol_pme2,
    names_to = "lhipp_method",
    values_to = "lhipp_volume"
  ) %>%
  pivot_longer(
    rhipp_data_vol2:rhipp_vol_pme2,
    names_to = "rhipp_method",
    values_to = "rhipp_volume"
  ) %>%
  mutate(
    lhipp_method = gsub("lhipp_", "", lhipp_method),
    rhipp_method = gsub("rhipp_", "", rhipp_method)
  ) %>%
  filter(lhipp_method == rhipp_method) %>%
  select(-rhipp_method) %>%
  separate(lhipp_method, c("source", "method"), sep = "_") %>%
  mutate(
    source_mod = ifelse(source == "data", source, method),
    method_mod = ifelse(source == "vol", source, method)
  ) %>%
  select(-source, -method) %>%
  rename(source = source_mod, method = method_mod)

hipp_long$method <- ifelse(grepl("pme1", hipp_long$source), "vol1", "vol2")
hipp_long$source <- gsub("pme1", "pme", hipp_long$source)
hipp_long$source <- gsub("pme2", "pme", hipp_long$source)

thal_long <- thal_info %>%
  pivot_longer(
    lthal_data_vol2:lthal_vol_pme2,
    names_to = "lthal_method",
    values_to = "lthal_volume"
  ) %>%
  pivot_longer(
    rthal_data_vol2:rthal_vol_pme2,
    names_to = "rthal_method",
    values_to = "rthal_volume"
  ) %>%
  mutate(
    lthal_method = gsub("lthal_", "", lthal_method),
    rthal_method = gsub("rthal_", "", rthal_method)
  ) %>%
  filter(lthal_method == rthal_method) %>%
  select(-rthal_method) %>%
  separate(lthal_method, c("source", "method"), sep = "_") %>%
  mutate(
    source_mod = ifelse(source == "data", source, method),
    method_mod = ifelse(source == "vol", source, method)
  ) %>%
  select(-source, -method) %>%
  rename(source = source_mod, method = method_mod)

thal_long$method <- ifelse(grepl("pme1", thal_long$source), "vol1", "vol2")
thal_long$source <- gsub("pme1", "pme", thal_long$source)
thal_long$source <- gsub("pme2", "pme", thal_long$source)

thal_info_ts <- as_tsibble(
  thal_info,
  key = patno,
  index = time_from_bl,
  regular = FALSE
)

hipp_info_ts <- as_tsibble(
  hipp_info,
  key = patno,
  index = time_from_bl,
  regular = FALSE
)


hipp_info_ts %>%
  sample_n_keys(3) %>%
  ggplot(aes(x = time_from_bl, group = patno, color = patno)) +
  geom_line(aes(y = lhipp_data_vol2), linetype = "solid") +
  geom_line(aes(y = lhipp_vol_lpme2), linetype = "dashed") +
  geom_line(aes(y = lhipp_vol_pme2), linetype = "dotted") +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Hippocampus Volume")

lhipp_data_vol_pred <- vector()
lhipp_lpme_vol_pred <- vector()
lhipp_pme_vol_pred <- vector()
rhipp_data_vol_pred <- vector()
rhipp_lpme_vol_pred <- vector()
rhipp_pme_vol_pred <- vector()

for (patno_val in unique(hipp_info$patno)) {
  patno_data <- hipp_info %>%
    filter(patno == patno_val)
  lhipp_data_lm <- lm(lhipp_data_vol2 ~ date, data = patno_data)
  lhipp_data_pred <- predict(lhipp_data_lm, newdata = patno_data)
  lhipp_data_vol_pred <- c(lhipp_data_vol_pred, lhipp_data_pred)

  lhipp_lpme_lm <- lm(lhipp_vol_lpme2 ~ date, data = patno_data)
  lhipp_lpme_pred <- predict(lhipp_lpme_lm, newdata = patno_data)
  lhipp_lpme_vol_pred <- c(lhipp_lpme_vol_pred, lhipp_lpme_pred)

  lhipp_pme_lm <- lm(lhipp_vol_pme2 ~ date, data = patno_data)
  lhipp_pme_pred <- predict(lhipp_pme_lm, newdata = patno_data)
  lhipp_pme_vol_pred <- c(lhipp_pme_vol_pred, lhipp_pme_pred)

  rhipp_data_lm <- lm(rhipp_data_vol2 ~ date, data = patno_data)
  rhipp_data_pred <- predict(rhipp_data_lm, newdata = patno_data)
  rhipp_data_vol_pred <- c(rhipp_data_vol_pred, rhipp_data_pred)

  rhipp_lpme_lm <- lm(rhipp_vol_lpme2 ~ date, data = patno_data)
  rhipp_lpme_pred <- predict(rhipp_lpme_lm, newdata = patno_data)
  rhipp_lpme_vol_pred <- c(rhipp_lpme_vol_pred, rhipp_lpme_pred)

  rhipp_pme_lm <- lm(rhipp_vol_pme2 ~ date, data = patno_data)
  rhipp_pme_pred <- predict(rhipp_pme_lm, newdata = patno_data)
  rhipp_pme_vol_pred <- c(rhipp_pme_vol_pred, rhipp_pme_pred)
}

hipp_lm_pred <- tibble(
  lhipp_data_vol_pred,
  lhipp_lpme_vol_pred,
  lhipp_pme_vol_pred,
  rhipp_data_vol_pred,
  rhipp_lpme_vol_pred,
  rhipp_pme_vol_pred
)

hipp_info <- bind_cols(hipp_info, hipp_lm_pred)

hipp_info_sd <- hipp_info %>%
  group_by(patno) %>%
  summarize(
    lhipp_data_sd = sd(lhipp_data_vol2),
    lhipp_lpme_sd = sd(lhipp_vol_lpme2),
    lhipp_pme_sd = sd(lhipp_vol_pme2),
    rhipp_data_sd = sd(rhipp_data_vol2),
    rhipp_lpme_sd = sd(rhipp_vol_lpme2),
    rhipp_pme_sd = sd(rhipp_vol_pme2),
    lhipp_data_adj_sd = sd(lhipp_data_vol2 - lhipp_data_vol_pred),
    lhipp_lpme_adj_sd = sd(lhipp_vol_lpme2 - lhipp_lpme_vol_pred),
    lhipp_pme_adj_sd = sd(lhipp_vol_pme2 - lhipp_pme_vol_pred),
    rhipp_data_adj_sd = sd(rhipp_data_vol2 - rhipp_data_vol_pred),
    rhipp_lpme_adj_sd = sd(rhipp_vol_lpme2 - rhipp_lpme_vol_pred),
    rhipp_pme_adj_sd = sd(rhipp_vol_pme2 - rhipp_pme_vol_pred)
  ) %>%
  ungroup()

hipp_sd_mean <- hipp_info_sd %>%
  summarize(
    lhipp_data_sd_mean = mean(lhipp_data_sd),
    lhipp_lpme_sd_mean = mean(lhipp_lpme_sd),
    lhipp_pme_sd_mean = mean(lhipp_pme_sd),
    rhipp_data_sd_mean = mean(rhipp_data_sd),
    rhipp_lpme_sd_mean = mean(rhipp_lpme_sd),
    rhipp_pme_sd_mean = mean(rhipp_pme_sd),
    lhipp_data_adj_sd_mean = mean(lhipp_data_adj_sd),
    lhipp_lpme_adj_sd_mean = mean(lhipp_lpme_adj_sd),
    lhipp_pme_adj_sd_mean = mean(lhipp_pme_adj_sd),
    rhipp_data_adj_sd_mean = mean(rhipp_data_adj_sd),
    rhipp_lpme_adj_sd_mean = mean(rhipp_lpme_adj_sd),
    rhipp_pme_adj_sd_mean = mean(rhipp_pme_adj_sd)
  )

print(hipp_sd_mean, width = Inf)

hipp_sd_med <- hipp_info_sd %>%
  summarize(
    lhipp_data_sd_med = median(lhipp_data_sd),
    lhipp_lpme_sd_med = median(lhipp_lpme_sd),
    lhipp_pme_sd_med = median(lhipp_pme_sd),
    rhipp_data_sd_med = median(rhipp_data_sd),
    rhipp_lpme_sd_med = median(rhipp_lpme_sd),
    rhipp_pme_sd_med = median(rhipp_pme_sd),
    lhipp_data_adj_sd_med = median(lhipp_data_adj_sd),
    lhipp_lpme_adj_sd_med = median(lhipp_lpme_adj_sd),
    lhipp_pme_adj_sd_med = median(lhipp_pme_adj_sd),
    rhipp_data_adj_sd_med = median(rhipp_data_adj_sd),
    rhipp_lpme_adj_sd_med = median(rhipp_lpme_adj_sd),
    rhipp_pme_adj_sd_med = median(rhipp_pme_adj_sd)
  )

print(hipp_sd_med, width = Inf)

thal_info_ts %>%
  sample_n_keys(3) %>%
  ggplot(aes(x = time_from_bl, group = patno, color = patno)) +
  geom_line(aes(y = lthal_data_vol2), linetype = "solid") +
  geom_line(aes(y = lthal_vol_lpme2), linetype = "dashed") +
  geom_line(aes(y = lthal_vol_pme2), linetype = "dotted") +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Thalamus Volume")

lthal_data_vol_pred <- vector()
lthal_lpme_vol_pred <- vector()
lthal_pme_vol_pred <- vector()
rthal_data_vol_pred <- vector()
rthal_lpme_vol_pred <- vector()
rthal_pme_vol_pred <- vector()

for (patno_val in unique(thal_info$patno)) {
  patno_data <- thal_info %>%
    filter(patno == patno_val)
  lthal_data_lm <- lm(lthal_data_vol2 ~ date, data = patno_data)
  lthal_data_pred <- predict(lthal_data_lm, newdata = patno_data)
  lthal_data_vol_pred <- c(lthal_data_vol_pred, lthal_data_pred)

  lthal_lpme_lm <- lm(lthal_vol_lpme2 ~ date, data = patno_data)
  lthal_lpme_pred <- predict(lthal_lpme_lm, newdata = patno_data)
  lthal_lpme_vol_pred <- c(lthal_lpme_vol_pred, lthal_lpme_pred)

  lthal_pme_lm <- lm(lthal_vol_pme2 ~ date, data = patno_data)
  lthal_pme_pred <- predict(lthal_pme_lm, newdata = patno_data)
  lthal_pme_vol_pred <- c(lthal_pme_vol_pred, lthal_pme_pred)

  rthal_data_lm <- lm(rthal_data_vol2 ~ date, data = patno_data)
  rthal_data_pred <- predict(rthal_data_lm, newdata = patno_data)
  rthal_data_vol_pred <- c(rthal_data_vol_pred, rthal_data_pred)

  rthal_lpme_lm <- lm(rthal_vol_lpme2 ~ date, data = patno_data)
  rthal_lpme_pred <- predict(rthal_lpme_lm, newdata = patno_data)
  rthal_lpme_vol_pred <- c(rthal_lpme_vol_pred, rthal_lpme_pred)

  rthal_pme_lm <- lm(rthal_vol_pme2 ~ date, data = patno_data)
  rthal_pme_pred <- predict(rthal_pme_lm, newdata = patno_data)
  rthal_pme_vol_pred <- c(rthal_pme_vol_pred, rthal_pme_pred)
}

thal_lm_pred <- tibble(
  lthal_data_vol_pred,
  lthal_lpme_vol_pred,
  lthal_pme_vol_pred,
  rthal_data_vol_pred,
  rthal_lpme_vol_pred,
  rthal_pme_vol_pred
)

thal_info <- bind_cols(thal_info, thal_lm_pred)


thal_info_sd <- thal_info %>%
  group_by(patno) %>%
  summarize(
    lthal_data_sd = sd(lthal_data_vol2),
    lthal_lpme_sd = sd(lthal_vol_lpme2),
    lthal_pme_sd = sd(lthal_vol_pme2),
    rthal_data_sd = sd(rthal_data_vol2),
    rthal_lpme_sd = sd(rthal_vol_lpme2),
    rthal_pme_sd = sd(rthal_vol_pme2),
    lthal_data_adj_sd = sd(lthal_data_vol2 - lthal_data_vol_pred),
    lthal_lpme_adj_sd = sd(lthal_vol_lpme2 - lthal_lpme_vol_pred),
    lthal_pme_adj_sd = sd(lthal_vol_pme2 - lthal_pme_vol_pred),
    rthal_data_adj_sd = sd(rthal_data_vol2 - rthal_data_vol_pred),
    rthal_lpme_adj_sd = sd(rthal_vol_lpme2 - rthal_lpme_vol_pred),
    rthal_pme_adj_sd = sd(rthal_vol_pme2 - rthal_pme_vol_pred)
  ) %>%
  ungroup()

thal_sd_mean <- thal_info_sd %>%
  summarize(
    lthal_data_sd_mean = mean(lthal_data_sd),
    lthal_lpme_sd_mean = mean(lthal_lpme_sd),
    lthal_pme_sd_mean = mean(lthal_pme_sd),
    rthal_data_sd_mean = mean(rthal_data_sd),
    rthal_lpme_sd_mean = mean(rthal_lpme_sd),
    rthal_pme_sd_mean = mean(rthal_pme_sd),
    lthal_data_adj_sd_mean = mean(lthal_data_adj_sd),
    lthal_lpme_adj_sd_mean = mean(lthal_lpme_adj_sd),
    lthal_pme_adj_sd_mean = mean(lthal_pme_adj_sd),
    rthal_data_adj_sd_mean = mean(rthal_data_adj_sd),
    rthal_lpme_adj_sd_mean = mean(rthal_lpme_adj_sd),
    rthal_pme_adj_sd_mean = mean(rthal_pme_adj_sd)
  )

print(thal_sd_mean, width = Inf)

thal_sd_med <- thal_info_sd %>%
  summarize(
    lthal_data_sd_med = median(lthal_data_sd),
    lthal_lpme_sd_med = median(lthal_lpme_sd),
    lthal_pme_sd_med = median(lthal_pme_sd),
    rthal_data_sd_med = median(rthal_data_sd),
    rthal_lpme_sd_med = median(rthal_lpme_sd),
    rthal_pme_sd_med = median(rthal_pme_sd),
    lthal_data_adj_sd_med = median(lthal_data_adj_sd),
    lthal_lpme_adj_sd_med = median(lthal_lpme_adj_sd),
    lthal_pme_adj_sd_med = median(lthal_pme_adj_sd),
    rthal_data_adj_sd_med = median(rthal_data_adj_sd),
    rthal_lpme_adj_sd_med = median(rthal_lpme_adj_sd),
    rthal_pme_adj_sd_med = median(rthal_pme_adj_sd)
  )

print(thal_sd_med, width = Inf)

set.seed(8310)
sampled_patnos <- hipp_info_ts %>%
  sample_n_keys(3) %>%
  .$patno %>%
  unique()

scale_factor <- 3
hipp_info_ts %>%
  filter(patno %in% sampled_patnos) %>%
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_point(aes(y = lhipp_data_vol2), shape = 16, color = colors[1]) +
  geom_line(aes(y = lhipp_data_vol2), color = colors[1]) +
  # geom_smooth(
  #   aes(y = lhipp_data_vol2),
  #   color = colors[1],
  #   method = "lm",
  #   se = FALSE,
  #   linetype = "dashed",
  #   linewidth = 1
  # ) +
  geom_point(aes(y = lhipp_vol_lpme2), shape = 15, color = colors[2]) +
  geom_line(aes(y = lhipp_vol_lpme2), color = colors[2]) +
  # geom_smooth(
  #   aes(y = lhipp_vol_lpme2),
  #   color = colors[2],
  #   method = "lm",
  #   se = FALSE,
  #   linetype = "dashed",
  #   linewidth = 1
  # ) +
  geom_point(aes(y = lhipp_vol_pme2), shape = 17, color = colors[3]) +
  geom_line(aes(y = lhipp_vol_pme2), color = colors[3]) +
  # geom_smooth(
  #   aes(y = lhipp_vol_pme2),
  #   color = colors[3],
  #   method = "lm",
  #   se = FALSE,
  #   linetype = "dashed",
  #   linewidth = 1
  # ) +
  facet_wrap(~patno) +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Hippocampus Volume") +
  theme(
    axis.title.x = element_text(size = scale_factor * 5),
    axis.title.y = element_text(size = scale_factor * 5),
    axis.text.x = element_text(size = scale_factor * 5),
    axis.text.y = element_text(size = scale_factor * 5),
    strip.text.x = element_text(size = scale_factor * 5)
  )
ggsave("paper/figures/adni_plots/adni_lhipp_volume_comp.png", dpi = 1500)

thal_info_ts %>%
  filter(patno %in% sampled_patnos) %>%
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_point(aes(y = lthal_data_vol2), shape = 16, color = colors[1]) +
  geom_line(aes(y = lthal_data_vol2), linetype = "solid", color = colors[1]) +
  # geom_smooth(
  #   aes(y = lthal_data_vol2),
  #   color = colors[1],
  #   method = "lm",
  #   se = FALSE,
  #   linetype = "dashed",
  #   linewidth = 1
  # ) +
  geom_point(aes(y = lthal_vol_lpme2), shape = 15, color = colors[2]) +
  geom_line(aes(y = lthal_vol_lpme2), linetype = "dashed", color = colors[2]) +
  # geom_smooth(
  #   aes(y = lthal_vol_lpme2),
  #   color = colors[2],
  #   method = "lm",
  #   se = FALSE,
  #   linetype = "dashed",
  #   linewidth = 1
  # ) +
  geom_point(aes(y = lthal_vol_pme2), shape = 17, color = colors[3]) +
  geom_line(aes(y = lthal_vol_pme2), linetype = "dotted", color = colors[3]) +
  # geom_smooth(
  #   aes(y = lthal_vol_pme2),
  #   color = colors[3],
  #   method = "lm",
  #   se = FALSE,
  #   linetype = "dashed",
  #   linewidth = 1
  # ) +
  facet_wrap(~patno) +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Thalamus Volume") +
  guides(
    color = guide_legend(title = "Participant ID"),
    shape = guide_legend(title = "Estimate Source")
  ) +
  theme(
    axis.title.x = element_text(size = scale_factor * 5),
    axis.title.y = element_text(size = scale_factor * 5),
    axis.text.x = element_text(size = scale_factor * 5),
    axis.text.y = element_text(size = scale_factor * 5),
    strip.text.x = element_text(size = scale_factor * 5)
  )
ggsave("paper/figures/adni_plots/adni_lthal_volume_comp.png", dpi = 1500)

adni_status %>%
  filter(
    PTID %in% eligible_patnos,
    VISCODE == "bl"
  ) %>%
  group_by(DX_bl) %>%
  tally()
