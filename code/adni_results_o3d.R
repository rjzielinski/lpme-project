library(dplyr)
library(here)
library(plotly)
library(pme)
library(reticulate)
library(RColorBrewer)

use_condaenv("lpme")
np <- import("numpy")
o3d <- import("open3d")
pv <- import("pyvista")

# source(here("code/functions/estimate_mesh_volume.R"))
source(here("code/functions/estimate_mesh_volume_poisson.R"))

result_dir <- here("output/adni")

result_files <- list.files(result_dir, full.names = TRUE, recursive = TRUE)

lhipp_files <- result_files[grepl("lhipp", result_files)]
rhipp_files <- result_files[grepl("rhipp", result_files)]
lthal_files <- result_files[grepl("lthal", result_files)]
rthal_files <- result_files[grepl("rthal", result_files)]

lthal_ex <- readRDS(lthal_files[1])

time_points <- unique(lthal_ex$data$time_from_bl)

palette <- colorRampPalette(c(brewer.pal(8, "Spectral")))
palette_colors <- palette(101)
time_colors <- palette_colors[ceiling((time_points + 1e-10) * 100)]

lthal_plot <- pv$Plotter(shape = tuple(4L, as.integer(length(time_points))))
lthal_sil_plot <- pv$Plotter(shape = tuple(2L, 2L))

for (time_idx in seq_along(time_points)) {
  temp_data <- lthal_ex$data |>
    filter(time_from_bl == time_points[time_idx]) |>
    select(x, y, z) |>
    as.matrix()

  temp_lpme_part_recon <- lthal_ex$lpme_part$reconstructions[
    lthal_ex$lpme_part$reconstructions[, 1] == time_points[time_idx],
    -1
  ]
  temp_pme_part_recon <- lthal_ex$pme_part$reconstructions[
    lthal_ex$pme_part$reconstructions[, 1] == time_points[time_idx],
    -1
  ]
  temp_pc_part_recon <- lthal_ex$pc_part$reconstructions[
    lthal_ex$pc_part$reconstructions[, 1] == time_points[time_idx],
    -1
  ]

  data_mesh <- estimate_mesh_volume_poisson(temp_data)
  lpme_mesh <- estimate_mesh_volume_poisson(temp_lpme_part_recon)
  pme_mesh <- estimate_mesh_volume_poisson(temp_pme_part_recon)
  pc_mesh <- estimate_mesh_volume_poisson(temp_pc_part_recon)

  lthal_plot$subplot(0L, as.integer(time_idx - 1))
  lthal_plot$add_mesh(
    # data_mesh$mesh$decimate(0.9, preserve_topology = TRUE),
    data_mesh$mesh,
    color = "gray",
    show_edges = TRUE,
    opacity = 0.5
  )

  lthal_plot$add_text(
    paste("Time =", round(time_points[time_idx], 2)),
    position = "upper_left",
    font_size = 12
  )
  lthal_plot$camera_position <- "xy"

  if (time_idx == 1) {
    lthal_plot$add_text(
      "Data",
      position = "lower_left",
      font_size = 12
    )
  }

  lthal_sil_plot$subplot(0L, 0L)
  lthal_sil_plot$add_silhouette(
    data_mesh$mesh,
    color = time_colors[time_idx],
    decimate = 0.9
  )

  if (time_idx == 1) {
    lthal_sil_plot$add_text(
      "Data",
      position = "lower_left",
      font_size = 12
    )
  }

  lthal_plot$subplot(1L, as.integer(time_idx - 1))
  lthal_plot$add_mesh(
    lpme_mesh$mesh,
    color = "blue",
    show_edges = TRUE,
    opacity = 0.75
  )
  lthal_plot$camera_position <- "xy"

  if (time_idx == 1) {
    lthal_plot$add_text(
      "LPME",
      position = "lower_left",
      font_size = 12
    )
  }

  lthal_sil_plot$subplot(0L, 1L)
  lthal_sil_plot$add_silhouette(
    lpme_mesh$mesh,
    color = time_colors[time_idx],
    decimate = 0.9
  )

  if (time_idx == 1) {
    lthal_sil_plot$add_text(
      "LPME",
      position = "lower_left",
      font_size = 12
    )
  }

  lthal_plot$subplot(2L, as.integer(time_idx - 1))
  lthal_plot$add_mesh(
    pme_mesh$mesh,
    color = "green",
    show_edges = TRUE,
    opacity = 0.75
  )
  lthal_plot$camera_position <- "xy"

  if (time_idx == 1) {
    lthal_plot$add_text(
      "PME",
      position = "lower_left",
      font_size = 12
    )
  }

  lthal_sil_plot$subplot(1L, 0L)
  lthal_sil_plot$add_silhouette(
    pme_mesh$mesh,
    color = time_colors[time_idx],
    decimate = 0.9
  )

  if (time_idx == 1) {
    lthal_sil_plot$add_text(
      "PME",
      position = "lower_left",
      font_size = 12
    )
  }

  lthal_plot$subplot(3L, as.integer(time_idx - 1))
  lthal_plot$add_mesh(
    pc_mesh$mesh,
    color = "purple",
    show_edges = TRUE,
    opacity = 0.75
  )
  lthal_plot$camera_position <- "xy"

  if (time_idx == 1) {
    lthal_plot$add_text(
      "PS",
      position = "lower_left",
      font_size = 12
    )
  }

  lthal_sil_plot$subplot(1L, 1L)
  lthal_sil_plot$add_silhouette(
    pc_mesh$mesh,
    color = time_colors[time_idx],
    decimate = 0.9
  )

  if (time_idx == 1) {
    lthal_sil_plot$add_text(
      "PS",
      position = "lower_left",
      font_size = 12
    )
  }
}

lthal_plot$add_title(
  "Example LPME Reconstruction of Left Thalamus",
  font_size = 24
)

lthal_plot$link_views()
lthal_plot$show()
lthal_plot$screenshot("test_lthal_plot.png")


lthal_sil_plot$link_views()
lthal_sil_plot$show()
lthal_sil_plot$screenshot("test_lthal_sil_plot.png")
