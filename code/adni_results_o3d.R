library(brolgar)
library(colorspace)
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

result_dir <- here("output/adni_area")

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

lhipp_idx <- which(grepl(common_patnos[17], lhipp_files))
lthal_idx <- which(grepl(common_patnos[17], lthal_files))

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

colors <- brewer.pal(4, "Set1")

adni_info <- read_csv(here("output/adni_areas.csv")) |>
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
    lhipp_data_dim1_sd = sd(lhipp_data_area_dim1, na.rm = TRUE),
    lhipp_data_dim2_sd = sd(lhipp_data_area_dim2, na.rm = TRUE),
    lhipp_data_dim3_sd = sd(lhipp_data_area_dim3, na.rm = TRUE),
    lhipp_lpme_dim1_sd = sd(lhipp_lpme_area_dim1, na.rm = TRUE),
    lhipp_lpme_dim2_sd = sd(lhipp_lpme_area_dim2, na.rm = TRUE),
    lhipp_lpme_dim3_sd = sd(lhipp_lpme_area_dim3, na.rm = TRUE),
    lhipp_pme_dim1_sd = sd(lhipp_pme_area_dim1, na.rm = TRUE),
    lhipp_pme_dim2_sd = sd(lhipp_pme_area_dim2, na.rm = TRUE),
    lhipp_pme_dim3_sd = sd(lhipp_pme_area_dim3, na.rm = TRUE),
    lhipp_pc_dim1_sd = sd(lhipp_pc_area_dim1, na.rm = TRUE),
    lhipp_pc_dim2_sd = sd(lhipp_pc_area_dim2, na.rm = TRUE),
    lhipp_pc_dim3_sd = sd(lhipp_pc_area_dim3, na.rm = TRUE),

    rhipp_data_dim1_sd = sd(rhipp_data_area_dim1, na.rm = TRUE),
    rhipp_data_dim2_sd = sd(rhipp_data_area_dim2, na.rm = TRUE),
    rhipp_data_dim3_sd = sd(rhipp_data_area_dim3, na.rm = TRUE),
    rhipp_lpme_dim1_sd = sd(rhipp_lpme_area_dim1, na.rm = TRUE),
    rhipp_lpme_dim2_sd = sd(rhipp_lpme_area_dim2, na.rm = TRUE),
    rhipp_lpme_dim3_sd = sd(rhipp_lpme_area_dim3, na.rm = TRUE),
    rhipp_pme_dim1_sd = sd(rhipp_pme_area_dim1, na.rm = TRUE),
    rhipp_pme_dim2_sd = sd(rhipp_pme_area_dim2, na.rm = TRUE),
    rhipp_pme_dim3_sd = sd(rhipp_pme_area_dim3, na.rm = TRUE),
    rhipp_pc_dim1_sd = sd(rhipp_pc_area_dim1, na.rm = TRUE),
    rhipp_pc_dim2_sd = sd(rhipp_pc_area_dim2, na.rm = TRUE),
    rhipp_pc_dim3_sd = sd(rhipp_pc_area_dim3, na.rm = TRUE)
  )


thal_summary <- adni_info %>%
  group_by(patno) %>%
  summarize(
    lthal_data_dim1_sd = sd(lthal_data_area_dim1, na.rm = TRUE),
    lthal_data_dim2_sd = sd(lthal_data_area_dim2, na.rm = TRUE),
    lthal_data_dim3_sd = sd(lthal_data_area_dim3, na.rm = TRUE),
    lthal_lpme_dim1_sd = sd(lthal_lpme_area_dim1, na.rm = TRUE),
    lthal_lpme_dim2_sd = sd(lthal_lpme_area_dim2, na.rm = TRUE),
    lthal_lpme_dim3_sd = sd(lthal_lpme_area_dim3, na.rm = TRUE),
    lthal_pme_dim1_sd = sd(lthal_pme_area_dim1, na.rm = TRUE),
    lthal_pme_dim2_sd = sd(lthal_pme_area_dim2, na.rm = TRUE),
    lthal_pme_dim3_sd = sd(lthal_pme_area_dim3, na.rm = TRUE),
    lthal_pc_dim1_sd = sd(lthal_pc_area_dim1, na.rm = TRUE),
    lthal_pc_dim2_sd = sd(lthal_pc_area_dim2, na.rm = TRUE),
    lthal_pc_dim3_sd = sd(lthal_pc_area_dim3, na.rm = TRUE),

    rthal_data_dim1_sd = sd(rthal_data_area_dim1, na.rm = TRUE),
    rthal_data_dim2_sd = sd(rthal_data_area_dim2, na.rm = TRUE),
    rthal_data_dim3_sd = sd(rthal_data_area_dim3, na.rm = TRUE),
    rthal_lpme_dim1_sd = sd(rthal_lpme_area_dim1, na.rm = TRUE),
    rthal_lpme_dim2_sd = sd(rthal_lpme_area_dim2, na.rm = TRUE),
    rthal_lpme_dim3_sd = sd(rthal_lpme_area_dim3, na.rm = TRUE),
    rthal_pme_dim1_sd = sd(rthal_pme_area_dim1, na.rm = TRUE),
    rthal_pme_dim2_sd = sd(rthal_pme_area_dim2, na.rm = TRUE),
    rthal_pme_dim3_sd = sd(rthal_pme_area_dim3, na.rm = TRUE),
    rthal_pc_dim1_sd = sd(rthal_pc_area_dim1, na.rm = TRUE),
    rthal_pc_dim2_sd = sd(rthal_pc_area_dim2, na.rm = TRUE),
    rthal_pc_dim3_sd = sd(rthal_pc_area_dim3, na.rm = TRUE)
  )


hipp_long <- adni_info %>%
  select(
    patno,
    date,
    lhipp_data_area_dim1:rhipp_pc_area_dim3
  ) |>
  pivot_longer(
    contains("lhipp") & contains("dim"),
    names_to = "lhipp_source",
    values_to = "lhipp_sil_area"
  ) %>%
  pivot_longer(
    contains("rhipp") & contains("dim"),
    names_to = "rhipp_source",
    values_to = "rhipp_sil_area"
  ) |>
  separate(
    lhipp_source,
    into = c(NA, "lhipp_source", NA, "lhipp_dim"),
    sep = "_"
  ) |>
  separate(
    rhipp_source,
    into = c(NA, "rhipp_source", NA, "rhipp_dim"),
    sep = "_"
  ) |>
  filter(lhipp_source == rhipp_source, lhipp_dim == rhipp_dim) %>%
  select(-rhipp_source, -rhipp_dim) |>
  rename(source = lhipp_source, dim = lhipp_dim) |>
  mutate(
    dim = gsub("dim", "", dim),
    dim = as.factor(dim)
  )

thal_long <- adni_info %>%
  select(
    patno,
    date,
    lthal_data_area_dim1:rthal_pc_area_dim3
  ) |>
  pivot_longer(
    contains("lthal") & contains("dim"),
    names_to = "lthal_source",
    values_to = "lthal_sil_area"
  ) %>%
  pivot_longer(
    contains("rthal") & contains("dim"),
    names_to = "rthal_source",
    values_to = "rthal_sil_area"
  ) %>%
  separate(
    lthal_source,
    into = c(NA, "lthal_source", NA, "lthal_dim"),
    sep = "_"
  ) |>
  separate(
    rthal_source,
    into = c(NA, "rthal_source", NA, "rthal_dim"),
    sep = "_"
  ) |>
  filter(lthal_source == rthal_source, lthal_dim == rthal_dim) %>%
  select(-rthal_source, -rthal_dim) |>
  rename(source = lthal_source, dim = lthal_dim) |>
  mutate(
    dim = gsub("dim", "", dim),
    dim = as.factor(dim)
  )

adni_info_ts <- as_tsibble(
  adni_info,
  key = patno,
  index = time_from_bl,
  regular = FALSE
)

adni_info_ts %>%
  sample_n_keys(3) %>%
  ggplot(aes(x = time_from_bl, group = patno, color = patno)) +
  geom_line(aes(y = lhipp_data_area_dim1), linetype = "solid") +
  geom_line(aes(y = lhipp_lpme_area_dim1), linetype = "dashed") +
  geom_line(aes(y = lhipp_pme_area_dim1), linetype = "dotted") +
  geom_line(aes(y = lhipp_pc_area_dim1), linetype = "dotdash") +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Hippocampus Silhouette Area (Dim 1)")

lhipp_data_pred_dim1 <- vector()
lhipp_data_pred_dim2 <- vector()
lhipp_data_pred_dim3 <- vector()
lhipp_lpme_pred_dim1 <- vector()
lhipp_lpme_pred_dim2 <- vector()
lhipp_lpme_pred_dim3 <- vector()
lhipp_pme_pred_dim1 <- vector()
lhipp_pme_pred_dim2 <- vector()
lhipp_pme_pred_dim3 <- vector()
lhipp_pc_pred_dim1 <- vector()
lhipp_pc_pred_dim2 <- vector()
lhipp_pc_pred_dim3 <- vector()

rhipp_data_pred_dim1 <- vector()
rhipp_data_pred_dim2 <- vector()
rhipp_data_pred_dim3 <- vector()
rhipp_lpme_pred_dim1 <- vector()
rhipp_lpme_pred_dim2 <- vector()
rhipp_lpme_pred_dim3 <- vector()
rhipp_pme_pred_dim1 <- vector()
rhipp_pme_pred_dim2 <- vector()
rhipp_pme_pred_dim3 <- vector()
rhipp_pc_pred_dim1 <- vector()
rhipp_pc_pred_dim2 <- vector()
rhipp_pc_pred_dim3 <- vector()

lthal_data_pred_dim1 <- vector()
lthal_data_pred_dim2 <- vector()
lthal_data_pred_dim3 <- vector()
lthal_lpme_pred_dim1 <- vector()
lthal_lpme_pred_dim2 <- vector()
lthal_lpme_pred_dim3 <- vector()
lthal_pme_pred_dim1 <- vector()
lthal_pme_pred_dim2 <- vector()
lthal_pme_pred_dim3 <- vector()
lthal_pc_pred_dim1 <- vector()
lthal_pc_pred_dim2 <- vector()
lthal_pc_pred_dim3 <- vector()

rthal_data_pred_dim1 <- vector()
rthal_data_pred_dim2 <- vector()
rthal_data_pred_dim3 <- vector()
rthal_lpme_pred_dim1 <- vector()
rthal_lpme_pred_dim2 <- vector()
rthal_lpme_pred_dim3 <- vector()
rthal_pme_pred_dim1 <- vector()
rthal_pme_pred_dim2 <- vector()
rthal_pme_pred_dim3 <- vector()
rthal_pc_pred_dim1 <- vector()
rthal_pc_pred_dim2 <- vector()
rthal_pc_pred_dim3 <- vector()

for (patno_val in unique(adni_info$patno)) {
  patno_data <- adni_info %>%
    filter(patno == patno_val)

  lhipp_data_lm_dim1 <- lm(lhipp_data_area_dim1 ~ date, data = patno_data)
  lhipp_data_pred1 <- predict(lhipp_data_lm_dim1, newdata = patno_data)
  lhipp_data_pred_dim1 <- c(lhipp_data_pred_dim1, lhipp_data_pred1)
  lhipp_data_lm_dim2 <- lm(lhipp_data_area_dim2 ~ date, data = patno_data)
  lhipp_data_pred2 <- predict(lhipp_data_lm_dim2, newdata = patno_data)
  lhipp_data_pred_dim2 <- c(lhipp_data_pred_dim2, lhipp_data_pred2)
  lhipp_data_lm_dim3 <- lm(lhipp_data_area_dim3 ~ date, data = patno_data)
  lhipp_data_pred3 <- predict(lhipp_data_lm_dim3, newdata = patno_data)
  lhipp_data_pred_dim3 <- c(lhipp_data_pred_dim3, lhipp_data_pred3)

  lhipp_lpme_lm_dim1 <- lm(lhipp_lpme_area_dim1 ~ date, data = patno_data)
  lhipp_lpme_pred1 <- predict(lhipp_lpme_lm_dim1, newdata = patno_data)
  lhipp_lpme_pred_dim1 <- c(lhipp_lpme_pred_dim1, lhipp_lpme_pred1)
  lhipp_lpme_lm_dim2 <- lm(lhipp_lpme_area_dim2 ~ date, data = patno_data)
  lhipp_lpme_pred2 <- predict(lhipp_lpme_lm_dim2, newdata = patno_data)
  lhipp_lpme_pred_dim2 <- c(lhipp_lpme_pred_dim2, lhipp_lpme_pred2)
  lhipp_lpme_lm_dim3 <- lm(lhipp_lpme_area_dim3 ~ date, data = patno_data)
  lhipp_lpme_pred3 <- predict(lhipp_lpme_lm_dim3, newdata = patno_data)
  lhipp_lpme_pred_dim3 <- c(lhipp_lpme_pred_dim3, lhipp_lpme_pred3)

  lhipp_pme_lm_dim1 <- lm(lhipp_pme_area_dim1 ~ date, data = patno_data)
  lhipp_pme_pred1 <- predict(lhipp_pme_lm_dim1, newdata = patno_data)
  lhipp_pme_pred_dim1 <- c(lhipp_pme_pred_dim1, lhipp_pme_pred1)
  lhipp_pme_lm_dim2 <- lm(lhipp_pme_area_dim2 ~ date, data = patno_data)
  lhipp_pme_pred2 <- predict(lhipp_pme_lm_dim2, newdata = patno_data)
  lhipp_pme_pred_dim2 <- c(lhipp_pme_pred_dim2, lhipp_pme_pred2)
  lhipp_pme_lm_dim3 <- lm(lhipp_pme_area_dim3 ~ date, data = patno_data)
  lhipp_pme_pred3 <- predict(lhipp_pme_lm_dim3, newdata = patno_data)
  lhipp_pme_pred_dim3 <- c(lhipp_pme_pred_dim3, lhipp_pme_pred3)

  lhipp_pc_lm_dim1 <- lm(lhipp_pc_area_dim1 ~ date, data = patno_data)
  lhipp_pc_pred1 <- predict(lhipp_pc_lm_dim1, newdata = patno_data)
  lhipp_pc_pred_dim1 <- c(lhipp_pc_pred_dim1, lhipp_pc_pred1)
  lhipp_pc_lm_dim2 <- lm(lhipp_pc_area_dim2 ~ date, data = patno_data)
  lhipp_pc_pred2 <- predict(lhipp_pc_lm_dim2, newdata = patno_data)
  lhipp_pc_pred_dim2 <- c(lhipp_pc_pred_dim2, lhipp_pc_pred2)
  lhipp_pc_lm_dim3 <- lm(lhipp_pc_area_dim3 ~ date, data = patno_data)
  lhipp_pc_pred3 <- predict(lhipp_pc_lm_dim3, newdata = patno_data)
  lhipp_pc_pred_dim3 <- c(lhipp_pc_pred_dim3, lhipp_pc_pred3)

  rhipp_data_lm_dim1 <- lm(rhipp_data_area_dim1 ~ date, data = patno_data)
  rhipp_data_pred1 <- predict(rhipp_data_lm_dim1, newdata = patno_data)
  rhipp_data_pred_dim1 <- c(rhipp_data_pred_dim1, rhipp_data_pred1)
  rhipp_data_lm_dim2 <- lm(rhipp_data_area_dim2 ~ date, data = patno_data)
  rhipp_data_pred2 <- predict(rhipp_data_lm_dim2, newdata = patno_data)
  rhipp_data_pred_dim2 <- c(rhipp_data_pred_dim2, rhipp_data_pred2)
  rhipp_data_lm_dim3 <- lm(rhipp_data_area_dim3 ~ date, data = patno_data)
  rhipp_data_pred3 <- predict(rhipp_data_lm_dim3, newdata = patno_data)
  rhipp_data_pred_dim3 <- c(rhipp_data_pred_dim3, rhipp_data_pred3)

  rhipp_lpme_lm_dim1 <- lm(rhipp_lpme_area_dim1 ~ date, data = patno_data)
  rhipp_lpme_pred1 <- predict(rhipp_lpme_lm_dim1, newdata = patno_data)
  rhipp_lpme_pred_dim1 <- c(rhipp_lpme_pred_dim1, rhipp_lpme_pred1)
  rhipp_lpme_lm_dim2 <- lm(rhipp_lpme_area_dim2 ~ date, data = patno_data)
  rhipp_lpme_pred2 <- predict(rhipp_lpme_lm_dim2, newdata = patno_data)
  rhipp_lpme_pred_dim2 <- c(rhipp_lpme_pred_dim2, rhipp_lpme_pred2)
  rhipp_lpme_lm_dim3 <- lm(rhipp_lpme_area_dim3 ~ date, data = patno_data)
  rhipp_lpme_pred3 <- predict(rhipp_lpme_lm_dim3, newdata = patno_data)
  rhipp_lpme_pred_dim3 <- c(rhipp_lpme_pred_dim3, rhipp_lpme_pred3)

  rhipp_pme_lm_dim1 <- lm(rhipp_pme_area_dim1 ~ date, data = patno_data)
  rhipp_pme_pred1 <- predict(rhipp_pme_lm_dim1, newdata = patno_data)
  rhipp_pme_pred_dim1 <- c(rhipp_pme_pred_dim1, rhipp_pme_pred1)
  rhipp_pme_lm_dim2 <- lm(rhipp_pme_area_dim2 ~ date, data = patno_data)
  rhipp_pme_pred2 <- predict(rhipp_pme_lm_dim2, newdata = patno_data)
  rhipp_pme_pred_dim2 <- c(rhipp_pme_pred_dim2, rhipp_pme_pred2)
  rhipp_pme_lm_dim3 <- lm(rhipp_pme_area_dim3 ~ date, data = patno_data)
  rhipp_pme_pred3 <- predict(rhipp_pme_lm_dim3, newdata = patno_data)
  rhipp_pme_pred_dim3 <- c(rhipp_pme_pred_dim3, rhipp_pme_pred3)

  rhipp_pc_lm_dim1 <- lm(rhipp_pc_area_dim1 ~ date, data = patno_data)
  rhipp_pc_pred1 <- predict(rhipp_pc_lm_dim1, newdata = patno_data)
  rhipp_pc_pred_dim1 <- c(rhipp_pc_pred_dim1, rhipp_pc_pred1)
  rhipp_pc_lm_dim2 <- lm(rhipp_pc_area_dim2 ~ date, data = patno_data)
  rhipp_pc_pred2 <- predict(rhipp_pc_lm_dim2, newdata = patno_data)
  rhipp_pc_pred_dim2 <- c(rhipp_pc_pred_dim2, rhipp_pc_pred2)
  rhipp_pc_lm_dim3 <- lm(rhipp_pc_area_dim3 ~ date, data = patno_data)
  rhipp_pc_pred3 <- predict(rhipp_pc_lm_dim3, newdata = patno_data)
  rhipp_pc_pred_dim3 <- c(rhipp_pc_pred_dim3, rhipp_pc_pred3)

  lthal_data_lm_dim1 <- lm(lthal_data_area_dim1 ~ date, data = patno_data)
  lthal_data_pred1 <- predict(lthal_data_lm_dim1, newdata = patno_data)
  lthal_data_pred_dim1 <- c(lthal_data_pred_dim1, lthal_data_pred1)
  lthal_data_lm_dim2 <- lm(lthal_data_area_dim2 ~ date, data = patno_data)
  lthal_data_pred2 <- predict(lthal_data_lm_dim2, newdata = patno_data)
  lthal_data_pred_dim2 <- c(lthal_data_pred_dim2, lthal_data_pred2)
  lthal_data_lm_dim3 <- lm(lthal_data_area_dim3 ~ date, data = patno_data)
  lthal_data_pred3 <- predict(lthal_data_lm_dim3, newdata = patno_data)
  lthal_data_pred_dim3 <- c(lthal_data_pred_dim3, lthal_data_pred3)

  lthal_lpme_lm_dim1 <- lm(lthal_lpme_area_dim1 ~ date, data = patno_data)
  lthal_lpme_pred1 <- predict(lthal_lpme_lm_dim1, newdata = patno_data)
  lthal_lpme_pred_dim1 <- c(lthal_lpme_pred_dim1, lthal_lpme_pred1)
  lthal_lpme_lm_dim2 <- lm(lthal_lpme_area_dim2 ~ date, data = patno_data)
  lthal_lpme_pred2 <- predict(lthal_lpme_lm_dim2, newdata = patno_data)
  lthal_lpme_pred_dim2 <- c(lthal_lpme_pred_dim2, lthal_lpme_pred2)
  lthal_lpme_lm_dim3 <- lm(lthal_lpme_area_dim3 ~ date, data = patno_data)
  lthal_lpme_pred3 <- predict(lthal_lpme_lm_dim3, newdata = patno_data)
  lthal_lpme_pred_dim3 <- c(lthal_lpme_pred_dim3, lthal_lpme_pred3)

  lthal_pme_lm_dim1 <- lm(lthal_pme_area_dim1 ~ date, data = patno_data)
  lthal_pme_pred1 <- predict(lthal_pme_lm_dim1, newdata = patno_data)
  lthal_pme_pred_dim1 <- c(lthal_pme_pred_dim1, lthal_pme_pred1)
  lthal_pme_lm_dim2 <- lm(lthal_pme_area_dim2 ~ date, data = patno_data)
  lthal_pme_pred2 <- predict(lthal_pme_lm_dim2, newdata = patno_data)
  lthal_pme_pred_dim2 <- c(lthal_pme_pred_dim2, lthal_pme_pred2)
  lthal_pme_lm_dim3 <- lm(lthal_pme_area_dim3 ~ date, data = patno_data)
  lthal_pme_pred3 <- predict(lthal_pme_lm_dim3, newdata = patno_data)
  lthal_pme_pred_dim3 <- c(lthal_pme_pred_dim3, lthal_pme_pred3)

  lthal_pc_lm_dim1 <- lm(lthal_pc_area_dim1 ~ date, data = patno_data)
  lthal_pc_pred1 <- predict(lthal_pc_lm_dim1, newdata = patno_data)
  lthal_pc_pred_dim1 <- c(lthal_pc_pred_dim1, lthal_pc_pred1)
  lthal_pc_lm_dim2 <- lm(lthal_pc_area_dim2 ~ date, data = patno_data)
  lthal_pc_pred2 <- predict(lthal_pc_lm_dim2, newdata = patno_data)
  lthal_pc_pred_dim2 <- c(lthal_pc_pred_dim2, lthal_pc_pred2)
  lthal_pc_lm_dim3 <- lm(lthal_pc_area_dim3 ~ date, data = patno_data)
  lthal_pc_pred3 <- predict(lthal_pc_lm_dim3, newdata = patno_data)
  lthal_pc_pred_dim3 <- c(lthal_pc_pred_dim3, lthal_pc_pred3)

  rthal_data_lm_dim1 <- lm(rthal_data_area_dim1 ~ date, data = patno_data)
  rthal_data_pred1 <- predict(rthal_data_lm_dim1, newdata = patno_data)
  rthal_data_pred_dim1 <- c(rthal_data_pred_dim1, rthal_data_pred1)
  rthal_data_lm_dim2 <- lm(rthal_data_area_dim2 ~ date, data = patno_data)
  rthal_data_pred2 <- predict(rthal_data_lm_dim2, newdata = patno_data)
  rthal_data_pred_dim2 <- c(rthal_data_pred_dim2, rthal_data_pred2)
  rthal_data_lm_dim3 <- lm(rthal_data_area_dim3 ~ date, data = patno_data)
  rthal_data_pred3 <- predict(rthal_data_lm_dim3, newdata = patno_data)
  rthal_data_pred_dim3 <- c(rthal_data_pred_dim3, rthal_data_pred3)

  rthal_lpme_lm_dim1 <- lm(rthal_lpme_area_dim1 ~ date, data = patno_data)
  rthal_lpme_pred1 <- predict(rthal_lpme_lm_dim1, newdata = patno_data)
  rthal_lpme_pred_dim1 <- c(rthal_lpme_pred_dim1, rthal_lpme_pred1)
  rthal_lpme_lm_dim2 <- lm(rthal_lpme_area_dim2 ~ date, data = patno_data)
  rthal_lpme_pred2 <- predict(rthal_lpme_lm_dim2, newdata = patno_data)
  rthal_lpme_pred_dim2 <- c(rthal_lpme_pred_dim2, rthal_lpme_pred2)
  rthal_lpme_lm_dim3 <- lm(rthal_lpme_area_dim3 ~ date, data = patno_data)
  rthal_lpme_pred3 <- predict(rthal_lpme_lm_dim3, newdata = patno_data)
  rthal_lpme_pred_dim3 <- c(rthal_lpme_pred_dim3, rthal_lpme_pred3)

  rthal_pme_lm_dim1 <- lm(rthal_pme_area_dim1 ~ date, data = patno_data)
  rthal_pme_pred1 <- predict(rthal_pme_lm_dim1, newdata = patno_data)
  rthal_pme_pred_dim1 <- c(rthal_pme_pred_dim1, rthal_pme_pred1)
  rthal_pme_lm_dim2 <- lm(rthal_pme_area_dim2 ~ date, data = patno_data)
  rthal_pme_pred2 <- predict(rthal_pme_lm_dim2, newdata = patno_data)
  rthal_pme_pred_dim2 <- c(rthal_pme_pred_dim2, rthal_pme_pred2)
  rthal_pme_lm_dim3 <- lm(rthal_pme_area_dim3 ~ date, data = patno_data)
  rthal_pme_pred3 <- predict(rthal_pme_lm_dim3, newdata = patno_data)
  rthal_pme_pred_dim3 <- c(rthal_pme_pred_dim3, rthal_pme_pred3)

  rthal_pc_lm_dim1 <- lm(rthal_pc_area_dim1 ~ date, data = patno_data)
  rthal_pc_pred1 <- predict(rthal_pc_lm_dim1, newdata = patno_data)
  rthal_pc_pred_dim1 <- c(rthal_pc_pred_dim1, rthal_pc_pred1)
  rthal_pc_lm_dim2 <- lm(rthal_pc_area_dim2 ~ date, data = patno_data)
  rthal_pc_pred2 <- predict(rthal_pc_lm_dim2, newdata = patno_data)
  rthal_pc_pred_dim2 <- c(rthal_pc_pred_dim2, rthal_pc_pred2)
  rthal_pc_lm_dim3 <- lm(rthal_pc_area_dim3 ~ date, data = patno_data)
  rthal_pc_pred3 <- predict(rthal_pc_lm_dim3, newdata = patno_data)
  rthal_pc_pred_dim3 <- c(rthal_pc_pred_dim3, rthal_pc_pred3)
}

adni_lm_pred <- tibble(
  lhipp_data_pred_dim1,
  lhipp_data_pred_dim2,
  lhipp_data_pred_dim3,
  lhipp_lpme_pred_dim1,
  lhipp_lpme_pred_dim2,
  lhipp_lpme_pred_dim3,
  lhipp_pme_pred_dim1,
  lhipp_pme_pred_dim2,
  lhipp_pme_pred_dim3,
  lhipp_pc_pred_dim1,
  lhipp_pc_pred_dim2,
  lhipp_pc_pred_dim3,

  rhipp_data_pred_dim1,
  rhipp_data_pred_dim2,
  rhipp_data_pred_dim3,
  rhipp_lpme_pred_dim1,
  rhipp_lpme_pred_dim2,
  rhipp_lpme_pred_dim3,
  rhipp_pme_pred_dim1,
  rhipp_pme_pred_dim2,
  rhipp_pme_pred_dim3,
  rhipp_pc_pred_dim1,
  rhipp_pc_pred_dim2,
  rhipp_pc_pred_dim3,

  lthal_data_pred_dim1,
  lthal_data_pred_dim2,
  lthal_data_pred_dim3,
  lthal_lpme_pred_dim1,
  lthal_lpme_pred_dim2,
  lthal_lpme_pred_dim3,
  lthal_pme_pred_dim1,
  lthal_pme_pred_dim2,
  lthal_pme_pred_dim3,
  lthal_pc_pred_dim1,
  lthal_pc_pred_dim2,
  lthal_pc_pred_dim3,

  rthal_data_pred_dim1,
  rthal_data_pred_dim2,
  rthal_data_pred_dim3,
  rthal_lpme_pred_dim1,
  rthal_lpme_pred_dim2,
  rthal_lpme_pred_dim3,
  rthal_pme_pred_dim1,
  rthal_pme_pred_dim2,
  rthal_pme_pred_dim3,
  rthal_pc_pred_dim1,
  rthal_pc_pred_dim2,
  rthal_pc_pred_dim3
)

adni_info <- bind_cols(adni_info, adni_lm_pred)

adni_info_sd <- adni_info %>%
  group_by(patno) %>%
  summarize(
    lhipp_data_sd_dim1 = sd(lhipp_data_area_dim1, na.rm = TRUE),
    lhipp_data_sd_dim2 = sd(lhipp_data_area_dim2, na.rm = TRUE),
    lhipp_data_sd_dim3 = sd(lhipp_data_area_dim3, na.rm = TRUE),
    lhipp_lpme_sd_dim1 = sd(lhipp_lpme_area_dim1, na.rm = TRUE),
    lhipp_lpme_sd_dim2 = sd(lhipp_lpme_area_dim2, na.rm = TRUE),
    lhipp_lpme_sd_dim3 = sd(lhipp_lpme_area_dim3, na.rm = TRUE),
    lhipp_pme_sd_dim1 = sd(lhipp_pme_area_dim1, na.rm = TRUE),
    lhipp_pme_sd_dim2 = sd(lhipp_pme_area_dim2, na.rm = TRUE),
    lhipp_pme_sd_dim3 = sd(lhipp_pme_area_dim3, na.rm = TRUE),
    lhipp_pc_sd_dim1 = sd(lhipp_pc_area_dim1, na.rm = TRUE),
    lhipp_pc_sd_dim2 = sd(lhipp_pc_area_dim2, na.rm = TRUE),
    lhipp_pc_sd_dim3 = sd(lhipp_pc_area_dim3, na.rm = TRUE),

    lhipp_data_adj_sd_dim1 = sd(
      lhipp_data_area_dim1 - lhipp_data_pred_dim1,
      na.rm = TRUE
    ),
    lhipp_data_adj_sd_dim2 = sd(
      lhipp_data_area_dim2 - lhipp_data_pred_dim2,
      na.rm = TRUE
    ),
    lhipp_data_adj_sd_dim3 = sd(
      lhipp_data_area_dim3 - lhipp_data_pred_dim3,
      na.rm = TRUE
    ),

    lhipp_lpme_adj_sd_dim1 = sd(
      lhipp_lpme_area_dim1 - lhipp_lpme_pred_dim1,
      na.rm = TRUE
    ),
    lhipp_lpme_adj_sd_dim2 = sd(
      lhipp_lpme_area_dim2 - lhipp_lpme_pred_dim2,
      na.rm = TRUE
    ),
    lhipp_lpme_adj_sd_dim3 = sd(
      lhipp_lpme_area_dim3 - lhipp_lpme_pred_dim3,
      na.rm = TRUE
    ),

    lhipp_pme_adj_sd_dim1 = sd(
      lhipp_pme_area_dim1 - lhipp_pme_pred_dim1,
      na.rm = TRUE
    ),
    lhipp_pme_adj_sd_dim2 = sd(
      lhipp_pme_area_dim2 - lhipp_pme_pred_dim2,
      na.rm = TRUE
    ),
    lhipp_pme_adj_sd_dim3 = sd(
      lhipp_pme_area_dim3 - lhipp_pme_pred_dim3,
      na.rm = TRUE
    ),

    lhipp_pc_adj_sd_dim1 = sd(
      lhipp_pc_area_dim1 - lhipp_pc_pred_dim1,
      na.rm = TRUE
    ),
    lhipp_pc_adj_sd_dim2 = sd(
      lhipp_pc_area_dim2 - lhipp_pc_pred_dim2,
      na.rm = TRUE
    ),
    lhipp_pc_adj_sd_dim3 = sd(
      lhipp_pc_area_dim3 - lhipp_pc_pred_dim3,
      na.rm = TRUE
    ),

    rhipp_data_sd_dim1 = sd(rhipp_data_area_dim1, na.rm = TRUE),
    rhipp_data_sd_dim2 = sd(rhipp_data_area_dim2, na.rm = TRUE),
    rhipp_data_sd_dim3 = sd(rhipp_data_area_dim3, na.rm = TRUE),
    rhipp_lpme_sd_dim1 = sd(rhipp_lpme_area_dim1, na.rm = TRUE),
    rhipp_lpme_sd_dim2 = sd(rhipp_lpme_area_dim2, na.rm = TRUE),
    rhipp_lpme_sd_dim3 = sd(rhipp_lpme_area_dim3, na.rm = TRUE),
    rhipp_pme_sd_dim1 = sd(rhipp_pme_area_dim1, na.rm = TRUE),
    rhipp_pme_sd_dim2 = sd(rhipp_pme_area_dim2, na.rm = TRUE),
    rhipp_pme_sd_dim3 = sd(rhipp_pme_area_dim3, na.rm = TRUE),
    rhipp_pc_sd_dim1 = sd(rhipp_pc_area_dim1, na.rm = TRUE),
    rhipp_pc_sd_dim2 = sd(rhipp_pc_area_dim2, na.rm = TRUE),
    rhipp_pc_sd_dim3 = sd(rhipp_pc_area_dim3, na.rm = TRUE),

    rhipp_data_adj_sd_dim1 = sd(
      rhipp_data_area_dim1 - rhipp_data_pred_dim1,
      na.rm = TRUE
    ),
    rhipp_data_adj_sd_dim2 = sd(
      rhipp_data_area_dim2 - rhipp_data_pred_dim2,
      na.rm = TRUE
    ),
    rhipp_data_adj_sd_dim3 = sd(
      rhipp_data_area_dim3 - rhipp_data_pred_dim3,
      na.rm = TRUE
    ),

    rhipp_lpme_adj_sd_dim1 = sd(
      rhipp_lpme_area_dim1 - rhipp_lpme_pred_dim1,
      na.rm = TRUE
    ),
    rhipp_lpme_adj_sd_dim2 = sd(
      rhipp_lpme_area_dim2 - rhipp_lpme_pred_dim2,
      na.rm = TRUE
    ),
    rhipp_lpme_adj_sd_dim3 = sd(
      rhipp_lpme_area_dim3 - rhipp_lpme_pred_dim3,
      na.rm = TRUE
    ),

    rhipp_pme_adj_sd_dim1 = sd(
      rhipp_pme_area_dim1 - rhipp_pme_pred_dim1,
      na.rm = TRUE
    ),
    rhipp_pme_adj_sd_dim2 = sd(
      rhipp_pme_area_dim2 - rhipp_pme_pred_dim2,
      na.rm = TRUE
    ),
    rhipp_pme_adj_sd_dim3 = sd(
      rhipp_pme_area_dim3 - rhipp_pme_pred_dim3,
      na.rm = TRUE
    ),

    rhipp_pc_adj_sd_dim1 = sd(
      rhipp_pc_area_dim1 - rhipp_pc_pred_dim1,
      na.rm = TRUE
    ),
    rhipp_pc_adj_sd_dim2 = sd(
      rhipp_pc_area_dim2 - rhipp_pc_pred_dim2,
      na.rm = TRUE
    ),
    rhipp_pc_adj_sd_dim3 = sd(
      rhipp_pc_area_dim3 - rhipp_pc_pred_dim3,
      na.rm = TRUE
    ),

    lthal_data_sd_dim1 = sd(lthal_data_area_dim1, na.rm = TRUE),
    lthal_data_sd_dim2 = sd(lthal_data_area_dim2, na.rm = TRUE),
    lthal_data_sd_dim3 = sd(lthal_data_area_dim3, na.rm = TRUE),
    lthal_lpme_sd_dim1 = sd(lthal_lpme_area_dim1, na.rm = TRUE),
    lthal_lpme_sd_dim2 = sd(lthal_lpme_area_dim2, na.rm = TRUE),
    lthal_lpme_sd_dim3 = sd(lthal_lpme_area_dim3, na.rm = TRUE),
    lthal_pme_sd_dim1 = sd(lthal_pme_area_dim1, na.rm = TRUE),
    lthal_pme_sd_dim2 = sd(lthal_pme_area_dim2, na.rm = TRUE),
    lthal_pme_sd_dim3 = sd(lthal_pme_area_dim3, na.rm = TRUE),
    lthal_pc_sd_dim1 = sd(lthal_pc_area_dim1, na.rm = TRUE),
    lthal_pc_sd_dim2 = sd(lthal_pc_area_dim2, na.rm = TRUE),
    lthal_pc_sd_dim3 = sd(lthal_pc_area_dim3, na.rm = TRUE),

    lthal_data_adj_sd_dim1 = sd(
      lthal_data_area_dim1 - lthal_data_pred_dim1,
      na.rm = TRUE
    ),
    lthal_data_adj_sd_dim2 = sd(
      lthal_data_area_dim2 - lthal_data_pred_dim2,
      na.rm = TRUE
    ),
    lthal_data_adj_sd_dim3 = sd(
      lthal_data_area_dim3 - lthal_data_pred_dim3,
      na.rm = TRUE
    ),

    lthal_lpme_adj_sd_dim1 = sd(
      lthal_lpme_area_dim1 - lthal_lpme_pred_dim1,
      na.rm = TRUE
    ),
    lthal_lpme_adj_sd_dim2 = sd(
      lthal_lpme_area_dim2 - lthal_lpme_pred_dim2,
      na.rm = TRUE
    ),
    lthal_lpme_adj_sd_dim3 = sd(
      lthal_lpme_area_dim3 - lthal_lpme_pred_dim3,
      na.rm = TRUE
    ),

    lthal_pme_adj_sd_dim1 = sd(
      lthal_pme_area_dim1 - lthal_pme_pred_dim1,
      na.rm = TRUE
    ),
    lthal_pme_adj_sd_dim2 = sd(
      lthal_pme_area_dim2 - lthal_pme_pred_dim2,
      na.rm = TRUE
    ),
    lthal_pme_adj_sd_dim3 = sd(
      lthal_pme_area_dim3 - lthal_pme_pred_dim3,
      na.rm = TRUE
    ),

    lthal_pc_adj_sd_dim1 = sd(
      lthal_pc_area_dim1 - lthal_pc_pred_dim1,
      na.rm = TRUE
    ),
    lthal_pc_adj_sd_dim2 = sd(
      lthal_pc_area_dim2 - lthal_pc_pred_dim2,
      na.rm = TRUE
    ),
    lthal_pc_adj_sd_dim3 = sd(
      lthal_pc_area_dim3 - lthal_pc_pred_dim3,
      na.rm = TRUE
    ),

    rthal_data_sd_dim1 = sd(rthal_data_area_dim1, na.rm = TRUE),
    rthal_data_sd_dim2 = sd(rthal_data_area_dim2, na.rm = TRUE),
    rthal_data_sd_dim3 = sd(rthal_data_area_dim3, na.rm = TRUE),
    rthal_lpme_sd_dim1 = sd(rthal_lpme_area_dim1, na.rm = TRUE),
    rthal_lpme_sd_dim2 = sd(rthal_lpme_area_dim2, na.rm = TRUE),
    rthal_lpme_sd_dim3 = sd(rthal_lpme_area_dim3, na.rm = TRUE),
    rthal_pme_sd_dim1 = sd(rthal_pme_area_dim1, na.rm = TRUE),
    rthal_pme_sd_dim2 = sd(rthal_pme_area_dim2, na.rm = TRUE),
    rthal_pme_sd_dim3 = sd(rthal_pme_area_dim3, na.rm = TRUE),
    rthal_pc_sd_dim1 = sd(rthal_pc_area_dim1, na.rm = TRUE),
    rthal_pc_sd_dim2 = sd(rthal_pc_area_dim2, na.rm = TRUE),
    rthal_pc_sd_dim3 = sd(rthal_pc_area_dim3, na.rm = TRUE),

    rthal_data_adj_sd_dim1 = sd(
      rthal_data_area_dim1 - rthal_data_pred_dim1,
      na.rm = TRUE
    ),
    rthal_data_adj_sd_dim2 = sd(
      rthal_data_area_dim2 - rthal_data_pred_dim2,
      na.rm = TRUE
    ),
    rthal_data_adj_sd_dim3 = sd(
      rthal_data_area_dim3 - rthal_data_pred_dim3,
      na.rm = TRUE
    ),

    rthal_lpme_adj_sd_dim1 = sd(
      rthal_lpme_area_dim1 - rthal_lpme_pred_dim1,
      na.rm = TRUE
    ),
    rthal_lpme_adj_sd_dim2 = sd(
      rthal_lpme_area_dim2 - rthal_lpme_pred_dim2,
      na.rm = TRUE
    ),
    rthal_lpme_adj_sd_dim3 = sd(
      rthal_lpme_area_dim3 - rthal_lpme_pred_dim3,
      na.rm = TRUE
    ),

    rthal_pme_adj_sd_dim1 = sd(
      rthal_pme_area_dim1 - rthal_pme_pred_dim1,
      na.rm = TRUE
    ),
    rthal_pme_adj_sd_dim2 = sd(
      rthal_pme_area_dim2 - rthal_pme_pred_dim2,
      na.rm = TRUE
    ),
    rthal_pme_adj_sd_dim3 = sd(
      rthal_pme_area_dim3 - rthal_pme_pred_dim3,
      na.rm = TRUE
    ),

    rthal_pc_adj_sd_dim1 = sd(
      rthal_pc_area_dim1 - rthal_pc_pred_dim1,
      na.rm = TRUE
    ),
    rthal_pc_adj_sd_dim2 = sd(
      rthal_pc_area_dim2 - rthal_pc_pred_dim2,
      na.rm = TRUE
    ),
    rthal_pc_adj_sd_dim3 = sd(
      rthal_pc_area_dim3 - rthal_pc_pred_dim3,
      na.rm = TRUE
    )
  ) %>%
  ungroup()

adni_sd_mean <- adni_info_sd %>%
  summarize(
    lhipp_data_sd_mean_dim1 = mean(lhipp_data_sd_dim1, na.rm = TRUE),
    lhipp_data_sd_mean_dim2 = mean(lhipp_data_sd_dim2, na.rm = TRUE),
    lhipp_data_sd_mean_dim3 = mean(lhipp_data_sd_dim3, na.rm = TRUE),
    lhipp_lpme_sd_mean_dim1 = mean(lhipp_lpme_sd_dim1, na.rm = TRUE),
    lhipp_lpme_sd_mean_dim2 = mean(lhipp_lpme_sd_dim2, na.rm = TRUE),
    lhipp_lpme_sd_mean_dim3 = mean(lhipp_lpme_sd_dim3, na.rm = TRUE),
    lhipp_pme_sd_mean_dim1 = mean(lhipp_pme_sd_dim1, na.rm = TRUE),
    lhipp_pme_sd_mean_dim2 = mean(lhipp_pme_sd_dim2, na.rm = TRUE),
    lhipp_pme_sd_mean_dim3 = mean(lhipp_pme_sd_dim3, na.rm = TRUE),
    lhipp_pc_sd_mean_dim1 = mean(lhipp_pc_sd_dim1, na.rm = TRUE),
    lhipp_pc_sd_mean_dim2 = mean(lhipp_pc_sd_dim2, na.rm = TRUE),
    lhipp_pc_sd_mean_dim3 = mean(lhipp_pc_sd_dim3, na.rm = TRUE),

    rhipp_data_sd_mean_dim1 = mean(rhipp_data_sd_dim1, na.rm = TRUE),
    rhipp_data_sd_mean_dim2 = mean(rhipp_data_sd_dim2, na.rm = TRUE),
    rhipp_data_sd_mean_dim3 = mean(rhipp_data_sd_dim3, na.rm = TRUE),
    rhipp_lpme_sd_mean_dim1 = mean(rhipp_lpme_sd_dim1, na.rm = TRUE),
    rhipp_lpme_sd_mean_dim2 = mean(rhipp_lpme_sd_dim2, na.rm = TRUE),
    rhipp_lpme_sd_mean_dim3 = mean(rhipp_lpme_sd_dim3, na.rm = TRUE),
    rhipp_pme_sd_mean_dim1 = mean(rhipp_pme_sd_dim1, na.rm = TRUE),
    rhipp_pme_sd_mean_dim2 = mean(rhipp_pme_sd_dim2, na.rm = TRUE),
    rhipp_pme_sd_mean_dim3 = mean(rhipp_pme_sd_dim3, na.rm = TRUE),
    rhipp_pc_sd_mean_dim1 = mean(rhipp_pc_sd_dim1, na.rm = TRUE),
    rhipp_pc_sd_mean_dim2 = mean(rhipp_pc_sd_dim2, na.rm = TRUE),
    rhipp_pc_sd_mean_dim3 = mean(rhipp_pc_sd_dim3, na.rm = TRUE),

    lhipp_data_adj_sd_mean_dim1 = mean(lhipp_data_adj_sd_dim1, na.rm = TRUE),
    lhipp_data_adj_sd_mean_dim2 = mean(lhipp_data_adj_sd_dim2, na.rm = TRUE),
    lhipp_data_adj_sd_mean_dim3 = mean(lhipp_data_adj_sd_dim3, na.rm = TRUE),
    lhipp_lpme_adj_sd_mean_dim1 = mean(lhipp_lpme_adj_sd_dim1, na.rm = TRUE),
    lhipp_lpme_adj_sd_mean_dim2 = mean(lhipp_lpme_adj_sd_dim2, na.rm = TRUE),
    lhipp_lpme_adj_sd_mean_dim3 = mean(lhipp_lpme_adj_sd_dim3, na.rm = TRUE),
    lhipp_pme_adj_sd_mean_dim1 = mean(lhipp_pme_adj_sd_dim1, na.rm = TRUE),
    lhipp_pme_adj_sd_mean_dim2 = mean(lhipp_pme_adj_sd_dim2, na.rm = TRUE),
    lhipp_pme_adj_sd_mean_dim3 = mean(lhipp_pme_adj_sd_dim3, na.rm = TRUE),
    lhipp_pc_adj_sd_mean_dim1 = mean(lhipp_pc_adj_sd_dim1, na.rm = TRUE),
    lhipp_pc_adj_sd_mean_dim2 = mean(lhipp_pc_adj_sd_dim2, na.rm = TRUE),
    lhipp_pc_adj_sd_mean_dim3 = mean(lhipp_pc_adj_sd_dim3, na.rm = TRUE),

    rhipp_data_adj_sd_mean_dim1 = mean(rhipp_data_adj_sd_dim1, na.rm = TRUE),
    rhipp_data_adj_sd_mean_dim2 = mean(rhipp_data_adj_sd_dim2, na.rm = TRUE),
    rhipp_data_adj_sd_mean_dim3 = mean(rhipp_data_adj_sd_dim3, na.rm = TRUE),
    rhipp_lpme_adj_sd_mean_dim1 = mean(rhipp_lpme_adj_sd_dim1, na.rm = TRUE),
    rhipp_lpme_adj_sd_mean_dim2 = mean(rhipp_lpme_adj_sd_dim2, na.rm = TRUE),
    rhipp_lpme_adj_sd_mean_dim3 = mean(rhipp_lpme_adj_sd_dim3, na.rm = TRUE),
    rhipp_pme_adj_sd_mean_dim1 = mean(rhipp_pme_adj_sd_dim1, na.rm = TRUE),
    rhipp_pme_adj_sd_mean_dim2 = mean(rhipp_pme_adj_sd_dim2, na.rm = TRUE),
    rhipp_pme_adj_sd_mean_dim3 = mean(rhipp_pme_adj_sd_dim3, na.rm = TRUE),
    rhipp_pc_adj_sd_mean_dim1 = mean(rhipp_pc_adj_sd_dim1, na.rm = TRUE),
    rhipp_pc_adj_sd_mean_dim2 = mean(rhipp_pc_adj_sd_dim2, na.rm = TRUE),
    rhipp_pc_adj_sd_mean_dim3 = mean(rhipp_pc_adj_sd_dim3, na.rm = TRUE),

    lthal_data_sd_mean_dim1 = mean(lthal_data_sd_dim1, na.rm = TRUE),
    lthal_data_sd_mean_dim2 = mean(lthal_data_sd_dim2, na.rm = TRUE),
    lthal_data_sd_mean_dim3 = mean(lthal_data_sd_dim3, na.rm = TRUE),
    lthal_lpme_sd_mean_dim1 = mean(lthal_lpme_sd_dim1, na.rm = TRUE),
    lthal_lpme_sd_mean_dim2 = mean(lthal_lpme_sd_dim2, na.rm = TRUE),
    lthal_lpme_sd_mean_dim3 = mean(lthal_lpme_sd_dim3, na.rm = TRUE),
    lthal_pme_sd_mean_dim1 = mean(lthal_pme_sd_dim1, na.rm = TRUE),
    lthal_pme_sd_mean_dim2 = mean(lthal_pme_sd_dim2, na.rm = TRUE),
    lthal_pme_sd_mean_dim3 = mean(lthal_pme_sd_dim3, na.rm = TRUE),
    lthal_pc_sd_mean_dim1 = mean(lthal_pc_sd_dim1, na.rm = TRUE),
    lthal_pc_sd_mean_dim2 = mean(lthal_pc_sd_dim2, na.rm = TRUE),
    lthal_pc_sd_mean_dim3 = mean(lthal_pc_sd_dim3, na.rm = TRUE),

    rthal_data_sd_mean_dim1 = mean(rthal_data_sd_dim1, na.rm = TRUE),
    rthal_data_sd_mean_dim2 = mean(rthal_data_sd_dim2, na.rm = TRUE),
    rthal_data_sd_mean_dim3 = mean(rthal_data_sd_dim3, na.rm = TRUE),
    rthal_lpme_sd_mean_dim1 = mean(rthal_lpme_sd_dim1, na.rm = TRUE),
    rthal_lpme_sd_mean_dim2 = mean(rthal_lpme_sd_dim2, na.rm = TRUE),
    rthal_lpme_sd_mean_dim3 = mean(rthal_lpme_sd_dim3, na.rm = TRUE),
    rthal_pme_sd_mean_dim1 = mean(rthal_pme_sd_dim1, na.rm = TRUE),
    rthal_pme_sd_mean_dim2 = mean(rthal_pme_sd_dim2, na.rm = TRUE),
    rthal_pme_sd_mean_dim3 = mean(rthal_pme_sd_dim3, na.rm = TRUE),
    rthal_pc_sd_mean_dim1 = mean(rthal_pc_sd_dim1, na.rm = TRUE),
    rthal_pc_sd_mean_dim2 = mean(rthal_pc_sd_dim2, na.rm = TRUE),
    rthal_pc_sd_mean_dim3 = mean(rthal_pc_sd_dim3, na.rm = TRUE),

    lthal_data_adj_sd_mean_dim1 = mean(lthal_data_adj_sd_dim1, na.rm = TRUE),
    lthal_data_adj_sd_mean_dim2 = mean(lthal_data_adj_sd_dim2, na.rm = TRUE),
    lthal_data_adj_sd_mean_dim3 = mean(lthal_data_adj_sd_dim3, na.rm = TRUE),
    lthal_lpme_adj_sd_mean_dim1 = mean(lthal_lpme_adj_sd_dim1, na.rm = TRUE),
    lthal_lpme_adj_sd_mean_dim2 = mean(lthal_lpme_adj_sd_dim2, na.rm = TRUE),
    lthal_lpme_adj_sd_mean_dim3 = mean(lthal_lpme_adj_sd_dim3, na.rm = TRUE),
    lthal_pme_adj_sd_mean_dim1 = mean(lthal_pme_adj_sd_dim1, na.rm = TRUE),
    lthal_pme_adj_sd_mean_dim2 = mean(lthal_pme_adj_sd_dim2, na.rm = TRUE),
    lthal_pme_adj_sd_mean_dim3 = mean(lthal_pme_adj_sd_dim3, na.rm = TRUE),
    lthal_pc_adj_sd_mean_dim1 = mean(lthal_pc_adj_sd_dim1, na.rm = TRUE),
    lthal_pc_adj_sd_mean_dim2 = mean(lthal_pc_adj_sd_dim2, na.rm = TRUE),
    lthal_pc_adj_sd_mean_dim3 = mean(lthal_pc_adj_sd_dim3, na.rm = TRUE),

    rthal_data_adj_sd_mean_dim1 = mean(rthal_data_adj_sd_dim1, na.rm = TRUE),
    rthal_data_adj_sd_mean_dim2 = mean(rthal_data_adj_sd_dim2, na.rm = TRUE),
    rthal_data_adj_sd_mean_dim3 = mean(rthal_data_adj_sd_dim3, na.rm = TRUE),
    rthal_lpme_adj_sd_mean_dim1 = mean(rthal_lpme_adj_sd_dim1, na.rm = TRUE),
    rthal_lpme_adj_sd_mean_dim2 = mean(rthal_lpme_adj_sd_dim2, na.rm = TRUE),
    rthal_lpme_adj_sd_mean_dim3 = mean(rthal_lpme_adj_sd_dim3, na.rm = TRUE),
    rthal_pme_adj_sd_mean_dim1 = mean(rthal_pme_adj_sd_dim1, na.rm = TRUE),
    rthal_pme_adj_sd_mean_dim2 = mean(rthal_pme_adj_sd_dim2, na.rm = TRUE),
    rthal_pme_adj_sd_mean_dim3 = mean(rthal_pme_adj_sd_dim3, na.rm = TRUE),
    rthal_pc_adj_sd_mean_dim1 = mean(rthal_pc_adj_sd_dim1, na.rm = TRUE),
    rthal_pc_adj_sd_mean_dim2 = mean(rthal_pc_adj_sd_dim2, na.rm = TRUE),
    rthal_pc_adj_sd_mean_dim3 = mean(rthal_pc_adj_sd_dim3, na.rm = TRUE),
  )

print(select(adni_sd_mean, contains("dim1")), width = Inf)
print(select(adni_sd_mean, contains("dim2")), width = Inf)
print(select(adni_sd_mean, contains("dim3")), width = Inf)


adni_sd_med <- adni_info_sd %>%
  summarize(
    lhipp_data_sd_median_dim1 = median(lhipp_data_sd_dim1, na.rm = TRUE),
    lhipp_data_sd_median_dim2 = median(lhipp_data_sd_dim2, na.rm = TRUE),
    lhipp_data_sd_median_dim3 = median(lhipp_data_sd_dim3, na.rm = TRUE),
    lhipp_lpme_sd_median_dim1 = median(lhipp_lpme_sd_dim1, na.rm = TRUE),
    lhipp_lpme_sd_median_dim2 = median(lhipp_lpme_sd_dim2, na.rm = TRUE),
    lhipp_lpme_sd_median_dim3 = median(lhipp_lpme_sd_dim3, na.rm = TRUE),
    lhipp_pme_sd_median_dim1 = median(lhipp_pme_sd_dim1, na.rm = TRUE),
    lhipp_pme_sd_median_dim2 = median(lhipp_pme_sd_dim2, na.rm = TRUE),
    lhipp_pme_sd_median_dim3 = median(lhipp_pme_sd_dim3, na.rm = TRUE),
    lhipp_pc_sd_median_dim1 = median(lhipp_pc_sd_dim1, na.rm = TRUE),
    lhipp_pc_sd_median_dim2 = median(lhipp_pc_sd_dim2, na.rm = TRUE),
    lhipp_pc_sd_median_dim3 = median(lhipp_pc_sd_dim3, na.rm = TRUE),

    rhipp_data_sd_median_dim1 = median(rhipp_data_sd_dim1, na.rm = TRUE),
    rhipp_data_sd_median_dim2 = median(rhipp_data_sd_dim2, na.rm = TRUE),
    rhipp_data_sd_median_dim3 = median(rhipp_data_sd_dim3, na.rm = TRUE),
    rhipp_lpme_sd_median_dim1 = median(rhipp_lpme_sd_dim1, na.rm = TRUE),
    rhipp_lpme_sd_median_dim2 = median(rhipp_lpme_sd_dim2, na.rm = TRUE),
    rhipp_lpme_sd_median_dim3 = median(rhipp_lpme_sd_dim3, na.rm = TRUE),
    rhipp_pme_sd_median_dim1 = median(rhipp_pme_sd_dim1, na.rm = TRUE),
    rhipp_pme_sd_median_dim2 = median(rhipp_pme_sd_dim2, na.rm = TRUE),
    rhipp_pme_sd_median_dim3 = median(rhipp_pme_sd_dim3, na.rm = TRUE),
    rhipp_pc_sd_median_dim1 = median(rhipp_pc_sd_dim1, na.rm = TRUE),
    rhipp_pc_sd_median_dim2 = median(rhipp_pc_sd_dim2, na.rm = TRUE),
    rhipp_pc_sd_median_dim3 = median(rhipp_pc_sd_dim3, na.rm = TRUE),

    lhipp_data_adj_sd_median_dim1 = median(
      lhipp_data_adj_sd_dim1,
      na.rm = TRUE
    ),
    lhipp_data_adj_sd_median_dim2 = median(
      lhipp_data_adj_sd_dim2,
      na.rm = TRUE
    ),
    lhipp_data_adj_sd_median_dim3 = median(
      lhipp_data_adj_sd_dim3,
      na.rm = TRUE
    ),
    lhipp_lpme_adj_sd_median_dim1 = median(
      lhipp_lpme_adj_sd_dim1,
      na.rm = TRUE
    ),
    lhipp_lpme_adj_sd_median_dim2 = median(
      lhipp_lpme_adj_sd_dim2,
      na.rm = TRUE
    ),
    lhipp_lpme_adj_sd_median_dim3 = median(
      lhipp_lpme_adj_sd_dim3,
      na.rm = TRUE
    ),

    lhipp_pme_adj_sd_median_dim1 = median(lhipp_pme_adj_sd_dim1, na.rm = TRUE),
    lhipp_pme_adj_sd_median_dim2 = median(lhipp_pme_adj_sd_dim2, na.rm = TRUE),
    lhipp_pme_adj_sd_median_dim3 = median(lhipp_pme_adj_sd_dim3, na.rm = TRUE),
    lhipp_pc_adj_sd_median_dim1 = median(lhipp_pc_adj_sd_dim1, na.rm = TRUE),
    lhipp_pc_adj_sd_median_dim2 = median(lhipp_pc_adj_sd_dim2, na.rm = TRUE),
    lhipp_pc_adj_sd_median_dim3 = median(lhipp_pc_adj_sd_dim3, na.rm = TRUE),

    rhipp_data_adj_sd_median_dim1 = median(
      rhipp_data_adj_sd_dim1,
      na.rm = TRUE
    ),
    rhipp_data_adj_sd_median_dim2 = median(
      rhipp_data_adj_sd_dim2,
      na.rm = TRUE
    ),
    rhipp_data_adj_sd_median_dim3 = median(
      rhipp_data_adj_sd_dim3,
      na.rm = TRUE
    ),
    rhipp_lpme_adj_sd_median_dim1 = median(
      rhipp_lpme_adj_sd_dim1,
      na.rm = TRUE
    ),
    rhipp_lpme_adj_sd_median_dim2 = median(
      rhipp_lpme_adj_sd_dim2,
      na.rm = TRUE
    ),
    rhipp_lpme_adj_sd_median_dim3 = median(
      rhipp_lpme_adj_sd_dim3,
      na.rm = TRUE
    ),

    rhipp_pme_adj_sd_median_dim1 = median(rhipp_pme_adj_sd_dim1, na.rm = TRUE),
    rhipp_pme_adj_sd_median_dim2 = median(rhipp_pme_adj_sd_dim2, na.rm = TRUE),
    rhipp_pme_adj_sd_median_dim3 = median(rhipp_pme_adj_sd_dim3, na.rm = TRUE),
    rhipp_pc_adj_sd_median_dim1 = median(rhipp_pc_adj_sd_dim1, na.rm = TRUE),
    rhipp_pc_adj_sd_median_dim2 = median(rhipp_pc_adj_sd_dim2, na.rm = TRUE),
    rhipp_pc_adj_sd_median_dim3 = median(rhipp_pc_adj_sd_dim3, na.rm = TRUE),

    lthal_data_sd_median_dim1 = median(lthal_data_sd_dim1, na.rm = TRUE),
    lthal_data_sd_median_dim2 = median(lthal_data_sd_dim2, na.rm = TRUE),
    lthal_data_sd_median_dim3 = median(lthal_data_sd_dim3, na.rm = TRUE),
    lthal_lpme_sd_median_dim1 = median(lthal_lpme_sd_dim1, na.rm = TRUE),
    lthal_lpme_sd_median_dim2 = median(lthal_lpme_sd_dim2, na.rm = TRUE),
    lthal_lpme_sd_median_dim3 = median(lthal_lpme_sd_dim3, na.rm = TRUE),
    lthal_pme_sd_median_dim1 = median(lthal_pme_sd_dim1, na.rm = TRUE),
    lthal_pme_sd_median_dim2 = median(lthal_pme_sd_dim2, na.rm = TRUE),
    lthal_pme_sd_median_dim3 = median(lthal_pme_sd_dim3, na.rm = TRUE),
    lthal_pc_sd_median_dim1 = median(lthal_pc_sd_dim1, na.rm = TRUE),
    lthal_pc_sd_median_dim2 = median(lthal_pc_sd_dim2, na.rm = TRUE),
    lthal_pc_sd_median_dim3 = median(lthal_pc_sd_dim3, na.rm = TRUE),

    rthal_data_sd_median_dim1 = median(rthal_data_sd_dim1, na.rm = TRUE),
    rthal_data_sd_median_dim2 = median(rthal_data_sd_dim2, na.rm = TRUE),
    rthal_data_sd_median_dim3 = median(rthal_data_sd_dim3, na.rm = TRUE),
    rthal_lpme_sd_median_dim1 = median(rthal_lpme_sd_dim1, na.rm = TRUE),
    rthal_lpme_sd_median_dim2 = median(rthal_lpme_sd_dim2, na.rm = TRUE),
    rthal_lpme_sd_median_dim3 = median(rthal_lpme_sd_dim3, na.rm = TRUE),
    rthal_pme_sd_median_dim1 = median(rthal_pme_sd_dim1, na.rm = TRUE),
    rthal_pme_sd_median_dim2 = median(rthal_pme_sd_dim2, na.rm = TRUE),
    rthal_pme_sd_median_dim3 = median(rthal_pme_sd_dim3, na.rm = TRUE),
    rthal_pc_sd_median_dim1 = median(rthal_pc_sd_dim1, na.rm = TRUE),
    rthal_pc_sd_median_dim2 = median(rthal_pc_sd_dim2, na.rm = TRUE),
    rthal_pc_sd_median_dim3 = median(rthal_pc_sd_dim3, na.rm = TRUE),

    lthal_data_adj_sd_median_dim1 = median(
      lthal_data_adj_sd_dim1,
      na.rm = TRUE
    ),
    lthal_data_adj_sd_median_dim2 = median(
      lthal_data_adj_sd_dim2,
      na.rm = TRUE
    ),
    lthal_data_adj_sd_median_dim3 = median(
      lthal_data_adj_sd_dim3,
      na.rm = TRUE
    ),
    lthal_lpme_adj_sd_median_dim1 = median(
      lthal_lpme_adj_sd_dim1,
      na.rm = TRUE
    ),
    lthal_lpme_adj_sd_median_dim2 = median(
      lthal_lpme_adj_sd_dim2,
      na.rm = TRUE
    ),
    lthal_lpme_adj_sd_median_dim3 = median(
      lthal_lpme_adj_sd_dim3,
      na.rm = TRUE
    ),

    lthal_pme_adj_sd_median_dim1 = median(lthal_pme_adj_sd_dim1, na.rm = TRUE),
    lthal_pme_adj_sd_median_dim2 = median(lthal_pme_adj_sd_dim2, na.rm = TRUE),
    lthal_pme_adj_sd_median_dim3 = median(lthal_pme_adj_sd_dim3, na.rm = TRUE),
    lthal_pc_adj_sd_median_dim1 = median(lthal_pc_adj_sd_dim1, na.rm = TRUE),
    lthal_pc_adj_sd_median_dim2 = median(lthal_pc_adj_sd_dim2, na.rm = TRUE),
    lthal_pc_adj_sd_median_dim3 = median(lthal_pc_adj_sd_dim3, na.rm = TRUE),

    rthal_data_adj_sd_median_dim1 = median(
      rthal_data_adj_sd_dim1,
      na.rm = TRUE
    ),
    rthal_data_adj_sd_median_dim2 = median(
      rthal_data_adj_sd_dim2,
      na.rm = TRUE
    ),
    rthal_data_adj_sd_median_dim3 = median(
      rthal_data_adj_sd_dim3,
      na.rm = TRUE
    ),
    rthal_lpme_adj_sd_median_dim1 = median(
      rthal_lpme_adj_sd_dim1,
      na.rm = TRUE
    ),
    rthal_lpme_adj_sd_median_dim2 = median(
      rthal_lpme_adj_sd_dim2,
      na.rm = TRUE
    ),
    rthal_lpme_adj_sd_median_dim3 = median(
      rthal_lpme_adj_sd_dim3,
      na.rm = TRUE
    ),

    rthal_pme_adj_sd_median_dim1 = median(rthal_pme_adj_sd_dim1, na.rm = TRUE),
    rthal_pme_adj_sd_median_dim2 = median(rthal_pme_adj_sd_dim2, na.rm = TRUE),
    rthal_pme_adj_sd_median_dim3 = median(rthal_pme_adj_sd_dim3, na.rm = TRUE),
    rthal_pc_adj_sd_median_dim1 = median(rthal_pc_adj_sd_dim1, na.rm = TRUE),
    rthal_pc_adj_sd_median_dim2 = median(rthal_pc_adj_sd_dim2, na.rm = TRUE),
    rthal_pc_adj_sd_median_dim3 = median(rthal_pc_adj_sd_dim3, na.rm = TRUE),
  )

print(select(adni_sd_med, contains("dim1")), width = Inf)
print(select(adni_sd_med, contains("dim2")), width = Inf)
print(select(adni_sd_med, contains("dim3")), width = Inf)

adni_info |>
  group_by(patno) |>
  summarize(
    lhipp_lpme_time = mean(lhipp_lpme_time),
    lhipp_pme_time = mean(lhipp_pme_time),
    lhipp_pc_time = mean(lhipp_pc_time),
    rhipp_lpme_time = mean(rhipp_lpme_time),
    rhipp_pme_time = mean(rhipp_pme_time),
    rhipp_pc_time = mean(rhipp_pc_time),
    lthal_lpme_time = mean(lthal_lpme_time),
    lthal_pme_time = mean(lthal_pme_time),
    lthal_pc_time = mean(lthal_pc_time),
    rthal_lpme_time = mean(rthal_lpme_time),
    rthal_pme_time = mean(rthal_pme_time),
    rthal_pc_time = mean(rthal_pc_time)
  ) |>
  ungroup() |>

  summarize(
    lhipp_lpme_mean = mean(lhipp_lpme_time, na.rm = TRUE),
    lhipp_pme_mean = mean(lhipp_pme_time, na.rm = TRUE),
    lhipp_pc_mean = mean(lhipp_pc_time, na.rm = TRUE),
    rhipp_lpme_mean = mean(rhipp_lpme_time, na.rm = TRUE),
    rhipp_pme_mean = mean(rhipp_pme_time, na.rm = TRUE),
    rhipp_pc_mean = mean(rhipp_pc_time, na.rm = TRUE),
    lthal_lpme_mean = mean(lthal_lpme_time, na.rm = TRUE),
    lthal_pme_mean = mean(lthal_pme_time, na.rm = TRUE),
    lthal_pc_mean = mean(lthal_pc_time, na.rm = TRUE),
    rthal_lpme_mean = mean(rthal_lpme_time, na.rm = TRUE),
    rthal_pme_mean = mean(rthal_pme_time, na.rm = TRUE),
    rthal_pc_mean = mean(rthal_pc_time, na.rm = TRUE)
  ) |>
  print(width = Inf)


set.seed(25)
sampled_patnos <- adni_info_ts %>%
  sample_n_keys(3) %>%
  .$patno %>%
  unique()

scale_factor <- 3
adni_info_ts %>%
  filter(patno %in% sampled_patnos) %>%
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_point(aes(y = lhipp_data_area_dim1), shape = 16, color = colors[1]) +
  geom_line(aes(y = lhipp_data_area_dim1), color = colors[1]) +
  geom_point(aes(y = lhipp_lpme_area_dim1), shape = 15, color = colors[2]) +
  geom_line(aes(y = lhipp_lpme_area_dim1), color = colors[2], linetype = 2) +
  geom_point(aes(y = lhipp_pme_area_dim1), shape = 17, color = colors[3]) +
  geom_line(aes(y = lhipp_pme_area_dim1), color = colors[3], linetype = 3) +
  geom_point(aes(y = lhipp_pc_area_dim1), shape = 18, color = colors[4]) +
  geom_line(aes(y = lhipp_pc_area_dim1), color = colors[4], linetype = 4) +

  facet_wrap(~patno) +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Hippocampus Profile Surface Area") +
  theme(
    axis.title.x = element_text(size = scale_factor * 5),
    axis.title.y = element_text(size = scale_factor * 5),
    axis.text.x = element_text(size = scale_factor * 5),
    axis.text.y = element_text(size = scale_factor * 5),
    strip.text.x = element_text(size = scale_factor * 5)
  )

ggsave("output/figures/adni_plots/adni_lhipp_area_comp_dim1.png", dpi = 1500)

scale_factor <- 3
adni_info_ts %>%
  filter(patno %in% sampled_patnos) %>%
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_point(aes(y = lhipp_data_area_dim2), shape = 16, color = colors[1]) +
  geom_line(aes(y = lhipp_data_area_dim2), color = colors[1]) +
  geom_point(aes(y = lhipp_lpme_area_dim2), shape = 15, color = colors[2]) +
  geom_line(aes(y = lhipp_lpme_area_dim2), color = colors[2], linetype = 2) +
  geom_point(aes(y = lhipp_pme_area_dim2), shape = 17, color = colors[3]) +
  geom_line(aes(y = lhipp_pme_area_dim2), color = colors[3], linetype = 3) +
  geom_point(aes(y = lhipp_pc_area_dim2), shape = 18, color = colors[4]) +
  geom_line(aes(y = lhipp_pc_area_dim2), color = colors[4], linetype = 4) +

  facet_wrap(~patno) +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Hippocampus Profile Surface Area") +
  theme(
    axis.title.x = element_text(size = scale_factor * 5),
    axis.title.y = element_text(size = scale_factor * 5),
    axis.text.x = element_text(size = scale_factor * 5),
    axis.text.y = element_text(size = scale_factor * 5),
    strip.text.x = element_text(size = scale_factor * 5)
  )

ggsave("output/figures/adni_plots/adni_lhipp_area_comp_dim2.png", dpi = 1500)

scale_factor <- 3
adni_info_ts %>%
  filter(patno %in% sampled_patnos) %>%
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_point(aes(y = lhipp_data_area_dim3), shape = 16, color = colors[1]) +
  geom_line(aes(y = lhipp_data_area_dim3), color = colors[1]) +
  geom_point(aes(y = lhipp_lpme_area_dim3), shape = 15, color = colors[2]) +
  geom_line(aes(y = lhipp_lpme_area_dim3), color = colors[2], linetype = 2) +
  geom_point(aes(y = lhipp_pme_area_dim3), shape = 17, color = colors[3]) +
  geom_line(aes(y = lhipp_pme_area_dim3), color = colors[3], linetype = 3) +
  geom_point(aes(y = lhipp_pc_area_dim3), shape = 18, color = colors[4]) +
  geom_line(aes(y = lhipp_pc_area_dim3), color = colors[4], linetype = 4) +

  facet_wrap(~patno) +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Hippocampus Profile Surface Area") +
  theme(
    axis.title.x = element_text(size = scale_factor * 5),
    axis.title.y = element_text(size = scale_factor * 5),
    axis.text.x = element_text(size = scale_factor * 5),
    axis.text.y = element_text(size = scale_factor * 5),
    strip.text.x = element_text(size = scale_factor * 5)
  )

ggsave("output/figures/adni_plots/adni_lhipp_area_comp_dim3.png", dpi = 1500)


scale_factor <- 3
adni_info_ts %>%
  filter(patno %in% sampled_patnos) %>%
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_point(aes(y = lthal_data_area_dim1), shape = 16, color = colors[1]) +
  geom_line(aes(y = lthal_data_area_dim1), color = colors[1]) +
  geom_point(aes(y = lthal_lpme_area_dim1), shape = 15, color = colors[2]) +
  geom_line(aes(y = lthal_lpme_area_dim1), color = colors[2], linetype = 2) +
  geom_point(aes(y = lthal_pme_area_dim1), shape = 17, color = colors[3]) +
  geom_line(aes(y = lthal_pme_area_dim1), color = colors[3], linetype = 3) +
  geom_point(aes(y = lthal_pc_area_dim1), shape = 18, color = colors[4]) +
  geom_line(aes(y = lthal_pc_area_dim1), color = colors[4], linetype = 4) +

  facet_wrap(~patno) +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Thalamus Profile Surface Area") +
  theme(
    axis.title.x = element_text(size = scale_factor * 5),
    axis.title.y = element_text(size = scale_factor * 5),
    axis.text.x = element_text(size = scale_factor * 5),
    axis.text.y = element_text(size = scale_factor * 5),
    strip.text.x = element_text(size = scale_factor * 5)
  )

ggsave("output/figures/adni_plots/adni_lthal_area_comp_dim1.png", dpi = 1500)


scale_factor <- 3
adni_info_ts %>%
  filter(patno %in% sampled_patnos) %>%
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_point(aes(y = lthal_data_area_dim2), shape = 16, color = colors[1]) +
  geom_line(aes(y = lthal_data_area_dim2), color = colors[1]) +
  geom_point(aes(y = lthal_lpme_area_dim2), shape = 15, color = colors[2]) +
  geom_line(aes(y = lthal_lpme_area_dim2), color = colors[2], linetype = 2) +
  geom_point(aes(y = lthal_pme_area_dim2), shape = 17, color = colors[3]) +
  geom_line(aes(y = lthal_pme_area_dim2), color = colors[3], linetype = 3) +
  geom_point(aes(y = lthal_pc_area_dim2), shape = 18, color = colors[4]) +
  geom_line(aes(y = lthal_pc_area_dim2), color = colors[4], linetype = 4) +

  facet_wrap(~patno) +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Thalamus Profile Surface Area") +
  theme(
    axis.title.x = element_text(size = scale_factor * 5),
    axis.title.y = element_text(size = scale_factor * 5),
    axis.text.x = element_text(size = scale_factor * 5),
    axis.text.y = element_text(size = scale_factor * 5),
    strip.text.x = element_text(size = scale_factor * 5)
  )

ggsave("output/figures/adni_plots/adni_lthal_area_comp_dim2.png", dpi = 1500)


scale_factor <- 3
adni_info_ts %>%
  filter(patno %in% sampled_patnos) %>%
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_point(aes(y = lthal_data_area_dim3), shape = 16, color = colors[1]) +
  geom_line(aes(y = lthal_data_area_dim3), color = colors[1]) +
  geom_point(aes(y = lthal_lpme_area_dim3), shape = 15, color = colors[2]) +
  geom_line(aes(y = lthal_lpme_area_dim3), color = colors[2], linetype = 2) +
  geom_point(aes(y = lthal_pme_area_dim3), shape = 17, color = colors[3]) +
  geom_line(aes(y = lthal_pme_area_dim3), color = colors[3], linetype = 3) +
  geom_point(aes(y = lthal_pc_area_dim3), shape = 18, color = colors[4]) +
  geom_line(aes(y = lthal_pc_area_dim3), color = colors[4], linetype = 4) +

  facet_wrap(~patno) +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Thalamus Profile Surface Area") +
  theme(
    axis.title.x = element_text(size = scale_factor * 5),
    axis.title.y = element_text(size = scale_factor * 5),
    axis.text.x = element_text(size = scale_factor * 5),
    axis.text.y = element_text(size = scale_factor * 5),
    strip.text.x = element_text(size = scale_factor * 5)
  )

ggsave("output/figures/adni_plots/adni_lthal_area_comp_dim3.png", dpi = 1500)

adni_status %>%
  filter(
    PTID %in% eligible_patnos,
    VISCODE == "bl"
  ) %>%
  group_by(DX_bl) %>%
  tally()
