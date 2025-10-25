library(brolgar)
library(dplyr)
library(here)
library(plotly)
library(pme)
library(readr)
library(reticulate)
library(RColorBrewer)
library(tidyverse)

use_condaenv("lpme")
np <- import("numpy")
o3d <- import("open3d")
pv <- import("pyvista")
# pymeshfix <- import("pymeshfix")

# source(here("code/functions/estimate_mesh_volume.R"))
source(here("code/functions/estimate_mesh_volume_poisson.R"))
source(here("code/functions/mesh_projection.R"))
source(here("code/functions/structure_plot_grid.R"))

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

common_patnos <- intersect(lhipp_patnos, lthal_patnos) |>
  intersect(rhipp_patnos) |>
  intersect(rthal_patnos)

lhipp_idx <- which(grepl(common_patnos[1], lhipp_files))
lthal_idx <- which(grepl(common_patnos[1], lthal_files))

lhipp_ex <- readRDS(lhipp_files[lhipp_idx])
lthal_ex <- readRDS(lthal_files[lthal_idx])

lhipp_plots <- structure_plot_grid(lhipp_ex, line_width = 5)

lhipp_plots$mesh_plot$show()
lhipp_plots$mesh_plot$screenshot(here("output/figures/lhipp_mesh_plot.png"))

lhipp_plots$sil_plots$dim1$show()
lhipp_plots$sil_plots$dim1$screenshot(here(
  "output/figures/lhipp_alg_silhouette_dim1.png"
))

lhipp_plots$sil_plots$dim2$show()
lhipp_plots$sil_plots$dim2$screenshot(here(
  "output/figures/lhipp_alg_silhouette_dim2.png"
))

lhipp_plots$sil_plots$dim3$show()
lhipp_plots$sil_plots$dim3$screenshot(here(
  "output/figures/lhipp_alg_silhouette_dim3.png"
))


lhipp_plots$time_sil_plots$dim1$show()
lhipp_plots$time_sil_plots$dim1$screenshot(here(
  "output/figures/lhipp_time_silhouette_dim1.png"
))

lhipp_plots$time_sil_plots$dim2$show()
lhipp_plots$time_sil_plots$dim2$screenshot(here(
  "output/figures/lhipp_time_silhouette_dim2.png"
))

lhipp_plots$time_sil_plots$dim3$show()
lhipp_plots$time_sil_plots$dim3$screenshot(here(
  "output/figures/lhipp_time_silhouette_dim3.png"
))


lthal_plots <- structure_plot_grid(
  lthal_ex,
  alpha_vals = exp(seq(-2, 0, 0.25)),
  line_width = 5
)

lthal_plots$mesh_plot$show()
lthal_plots$mesh_plot$screenshot(here("output/figures/lthal_mesh_plot.png"))

lthal_plots$sil_plots$dim1$show()
lthal_plots$sil_plots$dim1$screenshot(here(
  "output/figures/lthal_alg_silhouette_dim1.png"
))

lthal_plots$sil_plots$dim2$show()
lthal_plots$sil_plots$dim2$screenshot(here(
  "output/figures/lthal_alg_silhouette_dim2.png"
))

lthal_plots$sil_plots$dim3$show()
lthal_plots$sil_plots$dim3$screenshot(here(
  "output/figures/lthal_alg_silhouette_dim3.png"
))


lthal_plots$time_sil_plots$dim1$show()
lthal_plots$time_sil_plots$dim1$screenshot(here(
  "output/figures/lthal_time_silhouette_dim1.png"
))

lthal_plots$time_sil_plots$dim2$show()
lthal_plots$time_sil_plots$dim2$screenshot(here(
  "output/figures/lthal_time_silhouette_dim2.png"
))

lthal_plots$time_sil_plots$dim3$show()
lthal_plots$time_sil_plots$dim3$screenshot(here(
  "output/figures/lthal_time_silhouette_dim3.png"
))

colors <- brewer.pal(3, "Set1")

adni_info <- read_csv(here("output/adni_volumes.csv")) |>
  arrange(patno, date)

adni_bl <- adni_info |>
  group_by(patno) |>
  summarize(date_bl = min(date), max_date = max(date))

eligible_patnos <- adni_bl %>%
  filter(max_date - date_bl > 2) %>%
  select(patno) %>%
  unlist() %>%
  unique()

adni_info <- full_join(adni_info, adni_bl, by = "patno") |>
  mutate(time_from_bl = date - date_bl) |>
  filter(patno %in% eligible_patnos)

hipp_summary <- adni_info %>%
  group_by(patno) %>%
  summarize(
    lhipp_data_vol_sd = sd(lhipp_data_vol, na.rm = TRUE),
    lhipp_lpme_vol_sd = sd(lhipp_lpme_vol, na.rm = TRUE),
    lhipp_pme_vol_sd = sd(lhipp_pme_vol, na.rm = TRUE),
    lhipp_pc_vol_sd = sd(lhipp_pc_vol, na.rm = TRUE),
    rhipp_data_vol_sd = sd(rhipp_data_vol, na.rm = TRUE),
    rhipp_lpme_vol_sd = sd(rhipp_lpme_vol, na.rm = TRUE),
    rhipp_pme_vol_sd = sd(rhipp_pme_vol, na.rm = TRUE),
    rhipp_pc_vol_sd = sd(rhipp_pc_vol, na.rm = TRUE),
  )

thal_summary <- adni_info %>%
  group_by(patno) %>%
  summarize(
    lthal_data_vol_sd = sd(lthal_data_vol, na.rm = TRUE),
    lthal_lpme_vol_sd = sd(lthal_lpme_vol, na.rm = TRUE),
    lthal_pme_vol_sd = sd(lthal_pme_vol, na.rm = TRUE),
    lthal_pc_vol_sd = sd(lthal_pc_vol, na.rm = TRUE),
    rthal_data_vol_sd = sd(rthal_data_vol, na.rm = TRUE),
    rthal_lpme_vol_sd = sd(rthal_lpme_vol, na.rm = TRUE),
    rthal_pme_vol_sd = sd(rthal_pme_vol, na.rm = TRUE),
    rthal_pc_vol_sd = sd(rthal_pc_vol, na.rm = TRUE),
  )


hipp_long <- adni_info %>%
  select(
    patno,
    date,
    lhipp_data_vol:rhipp_pc_vol
  ) |>
  pivot_longer(
    lhipp_data_vol:lhipp_pc_vol,
    names_to = "lhipp_source",
    values_to = "lhipp_volume"
  ) %>%
  pivot_longer(
    rhipp_data_vol:rhipp_pc_vol,
    names_to = "rhipp_source",
    values_to = "rhipp_volume"
  ) %>%
  mutate(
    lhipp_source = gsub("lhipp_", "", lhipp_source),
    lhipp_source = gsub("_vol", "", lhipp_source),
    rhipp_source = gsub("rhipp_", "", rhipp_source),
    rhipp_source = gsub("_vol", "", rhipp_source)
  ) %>%
  filter(lhipp_source == rhipp_source) %>%
  select(-rhipp_source)

thal_long <- adni_info %>%
  select(
    patno,
    date,
    lthal_data_vol:rthal_pc_vol
  ) |>
  pivot_longer(
    lthal_data_vol:lthal_pc_vol,
    names_to = "lthal_source",
    values_to = "lthal_volume"
  ) %>%
  pivot_longer(
    rthal_data_vol:rthal_pc_vol,
    names_to = "rthal_source",
    values_to = "rthal_volume"
  ) %>%
  mutate(
    lthal_source = gsub("lthal_", "", lthal_source),
    lthal_source = gsub("_vol", "", lthal_source),
    rthal_source = gsub("rthal_", "", rthal_source),
    rthal_source = gsub("_vol", "", rthal_source)
  ) %>%
  filter(lthal_source == rthal_source) %>%
  select(-rthal_source)

adni_info_ts <- as_tsibble(
  adni_info,
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


adni_info_ts %>%
  sample_n_keys(3) %>%
  ggplot(aes(x = time_from_bl, group = patno, color = patno)) +
  geom_line(aes(y = lhipp_data_vol), linetype = "solid") +
  geom_line(aes(y = lhipp_lpme_vol), linetype = "dashed") +
  geom_line(aes(y = lhipp_pme_vol), linetype = "dotted") +
  geom_line(aes(y = lhipp_pc_vol), linetype = "dotdash") +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Hippocampus Volume")

lhipp_data_vol_pred <- vector()
lhipp_lpme_vol_pred <- vector()
lhipp_pme_vol_pred <- vector()
lhipp_pc_vol_pred <- vector()
rhipp_data_vol_pred <- vector()
rhipp_lpme_vol_pred <- vector()
rhipp_pme_vol_pred <- vector()
rhipp_pc_vol_pred <- vector()

lthal_data_vol_pred <- vector()
lthal_lpme_vol_pred <- vector()
lthal_pme_vol_pred <- vector()
lthal_pc_vol_pred <- vector()
rthal_data_vol_pred <- vector()
rthal_lpme_vol_pred <- vector()
rthal_pme_vol_pred <- vector()
rthal_pc_vol_pred <- vector()

for (patno_val in unique(adni_info$patno)) {
  patno_data <- adni_info %>%
    filter(patno == patno_val)
  lhipp_data_lm <- lm(lhipp_data_vol ~ date, data = patno_data)
  lhipp_data_pred <- predict(lhipp_data_lm, newdata = patno_data)
  lhipp_data_vol_pred <- c(lhipp_data_vol_pred, lhipp_data_pred)

  lhipp_lpme_lm <- lm(lhipp_lpme_vol ~ date, data = patno_data)
  lhipp_lpme_pred <- predict(lhipp_lpme_lm, newdata = patno_data)
  lhipp_lpme_vol_pred <- c(lhipp_lpme_vol_pred, lhipp_lpme_pred)

  lhipp_pme_lm <- lm(lhipp_pme_vol ~ date, data = patno_data)
  lhipp_pme_pred <- predict(lhipp_pme_lm, newdata = patno_data)
  lhipp_pme_vol_pred <- c(lhipp_pme_vol_pred, lhipp_pme_pred)

  lhipp_pc_lm <- lm(lhipp_pc_vol ~ date, data = patno_data)
  lhipp_pc_pred <- predict(lhipp_pc_lm, newdata = patno_data)
  lhipp_pc_vol_pred <- c(lhipp_pc_vol_pred, lhipp_pc_pred)

  rhipp_data_lm <- lm(rhipp_data_vol ~ date, data = patno_data)
  rhipp_data_pred <- predict(rhipp_data_lm, newdata = patno_data)
  rhipp_data_vol_pred <- c(rhipp_data_vol_pred, rhipp_data_pred)

  rhipp_lpme_lm <- lm(rhipp_lpme_vol ~ date, data = patno_data)
  rhipp_lpme_pred <- predict(rhipp_lpme_lm, newdata = patno_data)
  rhipp_lpme_vol_pred <- c(rhipp_lpme_vol_pred, rhipp_lpme_pred)

  rhipp_pme_lm <- lm(rhipp_pme_vol ~ date, data = patno_data)
  rhipp_pme_pred <- predict(rhipp_pme_lm, newdata = patno_data)
  rhipp_pme_vol_pred <- c(rhipp_pme_vol_pred, rhipp_pme_pred)

  rhipp_pc_lm <- lm(rhipp_pc_vol ~ date, data = patno_data)
  rhipp_pc_pred <- predict(rhipp_pc_lm, newdata = patno_data)
  rhipp_pc_vol_pred <- c(rhipp_pc_vol_pred, rhipp_pc_pred)

  lthal_data_lm <- lm(lthal_data_vol ~ date, data = patno_data)
  lthal_data_pred <- predict(lthal_data_lm, newdata = patno_data)
  lthal_data_vol_pred <- c(lthal_data_vol_pred, lthal_data_pred)

  lthal_lpme_lm <- lm(lthal_lpme_vol ~ date, data = patno_data)
  lthal_lpme_pred <- predict(lthal_lpme_lm, newdata = patno_data)
  lthal_lpme_vol_pred <- c(lthal_lpme_vol_pred, lthal_lpme_pred)

  lthal_pme_lm <- lm(lthal_pme_vol ~ date, data = patno_data)
  lthal_pme_pred <- predict(lthal_pme_lm, newdata = patno_data)
  lthal_pme_vol_pred <- c(lthal_pme_vol_pred, lthal_pme_pred)

  lthal_pc_lm <- lm(lthal_pc_vol ~ date, data = patno_data)
  lthal_pc_pred <- predict(lthal_pc_lm, newdata = patno_data)
  lthal_pc_vol_pred <- c(lthal_pc_vol_pred, lthal_pc_pred)

  rthal_data_lm <- lm(rthal_data_vol ~ date, data = patno_data)
  rthal_data_pred <- predict(rthal_data_lm, newdata = patno_data)
  rthal_data_vol_pred <- c(rthal_data_vol_pred, rthal_data_pred)

  rthal_lpme_lm <- lm(rthal_lpme_vol ~ date, data = patno_data)
  rthal_lpme_pred <- predict(rthal_lpme_lm, newdata = patno_data)
  rthal_lpme_vol_pred <- c(rthal_lpme_vol_pred, rthal_lpme_pred)

  rthal_pme_lm <- lm(rthal_pme_vol ~ date, data = patno_data)
  rthal_pme_pred <- predict(rthal_pme_lm, newdata = patno_data)
  rthal_pme_vol_pred <- c(rthal_pme_vol_pred, rthal_pme_pred)

  rthal_pc_lm <- lm(rthal_pc_vol ~ date, data = patno_data)
  rthal_pc_pred <- predict(rthal_pc_lm, newdata = patno_data)
  rthal_pc_vol_pred <- c(rthal_pc_vol_pred, rthal_pc_pred)
}

adni_lm_pred <- tibble(
  lhipp_data_vol_pred,
  lhipp_lpme_vol_pred,
  lhipp_pme_vol_pred,
  lhipp_pc_vol_pred,
  rhipp_data_vol_pred,
  rhipp_lpme_vol_pred,
  rhipp_pme_vol_pred,
  rhipp_pc_vol_pred,
  lthal_data_vol_pred,
  lthal_lpme_vol_pred,
  lthal_pme_vol_pred,
  lthal_pc_vol_pred,
  rthal_data_vol_pred,
  rthal_lpme_vol_pred,
  rthal_pme_vol_pred,
  rthal_pc_vol_pred
)

adni_info <- bind_cols(adni_info, adni_lm_pred)

adni_info_sd <- adni_info %>%
  group_by(patno) %>%
  summarize(
    lhipp_data_sd = sd(lhipp_data_vol, na.rm = TRUE),
    lhipp_lpme_sd = sd(lhipp_lpme_vol, na.rm = TRUE),
    lhipp_pme_sd = sd(lhipp_pme_vol, na.rm = TRUE),
    lhipp_pc_sd = sd(lhipp_pc_vol, na.rm = TRUE),

    rhipp_data_sd = sd(rhipp_data_vol, na.rm = TRUE),
    rhipp_lpme_sd = sd(rhipp_lpme_vol, na.rm = TRUE),
    rhipp_pme_sd = sd(rhipp_pme_vol, na.rm = TRUE),
    rhipp_pc_sd = sd(rhipp_pc_vol, na.rm = TRUE),

    lhipp_data_adj_sd = sd(lhipp_data_vol - lhipp_data_vol_pred, na.rm = TRUE),
    lhipp_lpme_adj_sd = sd(lhipp_lpme_vol - lhipp_lpme_vol_pred, na.rm = TRUE),
    lhipp_pme_adj_sd = sd(lhipp_pme_vol - lhipp_pme_vol_pred, na.rm = TRUE),
    lhipp_pc_adj_sd = sd(lhipp_pc_vol - lhipp_pc_vol_pred, na.rm = TRUE),

    rhipp_data_adj_sd = sd(rhipp_data_vol - rhipp_data_vol_pred, na.rm = TRUE),
    rhipp_lpme_adj_sd = sd(rhipp_lpme_vol - rhipp_lpme_vol_pred, na.rm = TRUE),
    rhipp_pme_adj_sd = sd(rhipp_pme_vol - rhipp_pme_vol_pred, na.rm = TRUE),
    rhipp_pc_adj_sd = sd(rhipp_pc_vol - rhipp_pc_vol_pred, na.rm = TRUE),

    lthal_data_sd = sd(lthal_data_vol, na.rm = TRUE),
    lthal_lpme_sd = sd(lthal_lpme_vol, na.rm = TRUE),
    lthal_pme_sd = sd(lthal_pme_vol, na.rm = TRUE),
    lthal_pc_sd = sd(lthal_pc_vol, na.rm = TRUE),

    rthal_data_sd = sd(rthal_data_vol, na.rm = TRUE),
    rthal_lpme_sd = sd(rthal_lpme_vol, na.rm = TRUE),
    rthal_pme_sd = sd(rthal_pme_vol, na.rm = TRUE),
    rthal_pc_sd = sd(rthal_pc_vol, na.rm = TRUE),

    lthal_data_adj_sd = sd(lthal_data_vol - lthal_data_vol_pred, na.rm = TRUE),
    lthal_lpme_adj_sd = sd(lthal_lpme_vol - lthal_lpme_vol_pred, na.rm = TRUE),
    lthal_pme_adj_sd = sd(lthal_pme_vol - lthal_pme_vol_pred, na.rm = TRUE),
    lthal_pc_adj_sd = sd(lthal_pc_vol - lthal_pc_vol_pred, na.rm = TRUE),

    rthal_data_adj_sd = sd(rthal_data_vol - rthal_data_vol_pred, na.rm = TRUE),
    rthal_lpme_adj_sd = sd(rthal_lpme_vol - rthal_lpme_vol_pred, na.rm = TRUE),
    rthal_pme_adj_sd = sd(rthal_pme_vol - rthal_pme_vol_pred, na.rm = TRUE),
    rthal_pc_adj_sd = sd(rthal_pc_vol - rthal_pc_vol_pred, na.rm = TRUE)
  ) %>%
  ungroup()

adni_sd_mean <- adni_info_sd %>%
  summarize(
    lhipp_data_sd_mean = mean(lhipp_data_sd, na.rm = TRUE),
    lhipp_lpme_sd_mean = mean(lhipp_lpme_sd, na.rm = TRUE),
    lhipp_pme_sd_mean = mean(lhipp_pme_sd, na.rm = TRUE),
    lhipp_pc_sd_mean = mean(lhipp_pc_sd, na.rm = TRUE),

    rhipp_data_sd_mean = mean(rhipp_data_sd, na.rm = TRUE),
    rhipp_lpme_sd_mean = mean(rhipp_lpme_sd, na.rm = TRUE),
    rhipp_pme_sd_mean = mean(rhipp_pme_sd, na.rm = TRUE),
    rhipp_pc_sd_mean = mean(rhipp_pc_sd, na.rm = TRUE),

    lhipp_data_adj_sd_mean = mean(lhipp_data_adj_sd, na.rm = TRUE),
    lhipp_lpme_adj_sd_mean = mean(lhipp_lpme_adj_sd, na.rm = TRUE),
    lhipp_pme_adj_sd_mean = mean(lhipp_pme_adj_sd, na.rm = TRUE),
    lhipp_pc_adj_sd_mean = mean(lhipp_pc_adj_sd, na.rm = TRUE),

    rhipp_data_adj_sd_mean = mean(rhipp_data_adj_sd, na.rm = TRUE),
    rhipp_lpme_adj_sd_mean = mean(rhipp_lpme_adj_sd, na.rm = TRUE),
    rhipp_pme_adj_sd_mean = mean(rhipp_pme_adj_sd, na.rm = TRUE),
    rhipp_pc_adj_sd_mean = mean(rhipp_pc_adj_sd, na.rm = TRUE),

    lthal_data_sd_mean = mean(lthal_data_sd, na.rm = TRUE),
    lthal_lpme_sd_mean = mean(lthal_lpme_sd, na.rm = TRUE),
    lthal_pme_sd_mean = mean(lthal_pme_sd, na.rm = TRUE),
    lthal_pc_sd_mean = mean(lthal_pc_sd, na.rm = TRUE),

    rthal_data_sd_mean = mean(rthal_data_sd, na.rm = TRUE),
    rthal_lpme_sd_mean = mean(rthal_lpme_sd, na.rm = TRUE),
    rthal_pme_sd_mean = mean(rthal_pme_sd, na.rm = TRUE),
    rthal_pc_sd_mean = mean(rthal_pc_sd, na.rm = TRUE),

    lthal_data_adj_sd_mean = mean(lthal_data_adj_sd, na.rm = TRUE),
    lthal_lpme_adj_sd_mean = mean(lthal_lpme_adj_sd, na.rm = TRUE),
    lthal_pme_adj_sd_mean = mean(lthal_pme_adj_sd, na.rm = TRUE),
    lthal_pc_adj_sd_mean = mean(lthal_pc_adj_sd, na.rm = TRUE),

    rthal_data_adj_sd_mean = mean(rthal_data_adj_sd, na.rm = TRUE),
    rthal_lpme_adj_sd_mean = mean(rthal_lpme_adj_sd, na.rm = TRUE),
    rthal_pme_adj_sd_mean = mean(rthal_pme_adj_sd, na.rm = TRUE),
    rthal_pc_adj_sd_mean = mean(rthal_pc_adj_sd, na.rm = TRUE)
  )

print(adni_sd_mean, width = Inf)

adni_sd_med <- adni_info_sd %>%
  summarize(
    lhipp_data_sd_med = median(lhipp_data_sd, na.rm = TRUE),
    lhipp_lpme_sd_med = median(lhipp_lpme_sd, na.rm = TRUE),
    lhipp_pme_sd_med = median(lhipp_pme_sd, na.rm = TRUE),
    lhipp_pc_sd_med = median(lhipp_pc_sd, na.rm = TRUE),

    rhipp_data_sd_med = median(rhipp_data_sd, na.rm = TRUE),
    rhipp_lpme_sd_med = median(rhipp_lpme_sd, na.rm = TRUE),
    rhipp_pme_sd_med = median(rhipp_pme_sd, na.rm = TRUE),
    rhipp_pc_sd_med = median(rhipp_pc_sd, na.rm = TRUE),

    lhipp_data_adj_sd_med = median(lhipp_data_adj_sd, na.rm = TRUE),
    lhipp_lpme_adj_sd_med = median(lhipp_lpme_adj_sd, na.rm = TRUE),
    lhipp_pme_adj_sd_med = median(lhipp_pme_adj_sd, na.rm = TRUE),
    lhipp_pc_adj_sd_med = median(lhipp_pc_adj_sd, na.rm = TRUE),

    rhipp_data_adj_sd_med = median(rhipp_data_adj_sd, na.rm = TRUE),
    rhipp_lpme_adj_sd_med = median(rhipp_lpme_adj_sd, na.rm = TRUE),
    rhipp_pme_adj_sd_med = median(rhipp_pme_adj_sd, na.rm = TRUE),
    rhipp_pc_adj_sd_med = median(rhipp_pc_adj_sd, na.rm = TRUE),

    lthal_data_sd_med = median(lthal_data_sd, na.rm = TRUE),
    lthal_lpme_sd_med = median(lthal_lpme_sd, na.rm = TRUE),
    lthal_pme_sd_med = median(lthal_pme_sd, na.rm = TRUE),
    lthal_pc_sd_med = median(lthal_pc_sd, na.rm = TRUE),

    rthal_data_sd_med = median(rthal_data_sd, na.rm = TRUE),
    rthal_lpme_sd_med = median(rthal_lpme_sd, na.rm = TRUE),
    rthal_pme_sd_med = median(rthal_pme_sd, na.rm = TRUE),
    rthal_pc_sd_med = median(rthal_pc_sd, na.rm = TRUE),

    lthal_data_adj_sd_med = median(lthal_data_adj_sd, na.rm = TRUE),
    lthal_lpme_adj_sd_med = median(lthal_lpme_adj_sd, na.rm = TRUE),
    lthal_pme_adj_sd_med = median(lthal_pme_adj_sd, na.rm = TRUE),
    lthal_pc_adj_sd_med = median(lthal_pc_adj_sd, na.rm = TRUE),

    rthal_data_adj_sd_med = median(rthal_data_adj_sd, na.rm = TRUE),
    rthal_lpme_adj_sd_med = median(rthal_lpme_adj_sd, na.rm = TRUE),
    rthal_pme_adj_sd_med = median(rthal_pme_adj_sd, na.rm = TRUE),
    rthal_pc_adj_sd_med = median(rthal_pc_adj_sd, na.rm = TRUE)
  )


print(adni_sd_med, width = Inf)

set.seed(412)
sampled_patnos <- adni_info_ts %>%
  sample_n_keys(3) %>%
  .$patno %>%
  unique()

scale_factor <- 3
adni_info_ts %>%
  filter(patno %in% sampled_patnos) %>%
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_point(aes(y = lhipp_data_vol), shape = 16, color = colors[1]) +
  geom_line(aes(y = lhipp_data_vol), color = colors[1]) +
  geom_point(aes(y = lhipp_lpme_vol), shape = 15, color = colors[2]) +
  geom_line(aes(y = lhipp_lpme_vol), color = colors[2]) +
  geom_point(aes(y = lhipp_pme_vol), shape = 17, color = colors[3]) +
  geom_line(aes(y = lhipp_pme_vol), color = colors[3]) +
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

scale_factor <- 3
adni_info_ts %>%
  filter(patno %in% sampled_patnos) %>%
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_point(aes(y = lthal_data_vol), shape = 16, color = colors[1]) +
  geom_line(aes(y = lthal_data_vol), color = colors[1]) +
  geom_point(aes(y = lthal_lpme_vol), shape = 15, color = colors[2]) +
  geom_line(aes(y = lthal_lpme_vol), color = colors[2]) +
  geom_point(aes(y = lthal_pme_vol), shape = 17, color = colors[3]) +
  geom_line(aes(y = lthal_pme_vol), color = colors[3]) +
  facet_wrap(~patno) +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left thalocampus Volume") +
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
