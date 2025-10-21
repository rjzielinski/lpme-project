structure_plot_grid <- function(
  adni_est,
  opacity = 0.75,
  line_width = 2
) {
  require("reticulate")
  use_condaenv("lpme")
  pv <- import("pyvista")

  require(here)

  source(here("code/functions/estimate_mesh_volume_poisson.R"))
  source(here("code/functions/mesh_projection.R"))
  time_points <- unique(adni_est$data$time_from_bl)

  palette <- colorRampPalette(c(brewer.pal(8, "Spectral")))
  palette_colors <- palette(101)
  time_colors <- palette_colors[ceiling((time_points + 1e-10) * 100)]

  adni_plot <- pv$Plotter(shape = tuple(4L, as.integer(length(time_points))))
  adni_sil_plot <- pv$Plotter(shape = tuple(2L, 2L))
  time_sil_plot <- pv$Plotter(
    shape = tuple(2L, as.integer(ceiling(length(time_points) / 2)))
  )

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

    data_projection <- mesh_projection(temp_data)
    lpme_projection <- mesh_projection(temp_lpme_part_recon)
    pme_projection <- mesh_projection(temp_pme_part_recon)
    pc_projection <- mesh_projection(temp_pc_part_recon)

    adni_plot$subplot(0L, as.integer(time_idx - 1))
    adni_plot$add_mesh(
      # data_mesh$mesh$decimate(0.9, preserve_topology = TRUE),
      data_mesh$mesh,
      color = "red",
      show_edges = TRUE,
      opacity = 0.5,
      label = "Data"
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
    adni_sil_plot$add_mesh(
      data_projection$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste0("Time = ", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot$add_text(
        "Data",
        position = "upper_left",
        font_size = 12
      )
    }
    if (time_idx == length(time_points)) {
      adni_sil_plot$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    time_sil_plot$subplot(
      as.integer(round(time_idx / length(time_points), 0)),
      as.integer((time_idx - 1) %% ceiling(length(time_points) / 2))
    )
    time_sil_plot$add_mesh(
      data_projection$boundary,
      show_edges = TRUE,
      color = "red",
      line_width = line_width,
      opacity = opacity,
      label = "Data"
    )
    time_sil_plot$add_text(
      paste("Time =", round(time_points[time_idx], 2)),
      position = "upper_left",
      font_size = 12
    )

    adni_plot$subplot(1L, as.integer(time_idx - 1))
    adni_plot$add_mesh(
      lpme_mesh$mesh,
      color = "blue",
      show_edges = TRUE,
      opacity = 0.75,
      label = "LPME"
    )

    if (time_idx == 1) {
      adni_plot$add_text(
        "LPME",
        position = "lower_left",
        font_size = 12
      )
    }

    adni_sil_plot$subplot(0L, 1L)
    adni_sil_plot$add_mesh(
      lpme_projection$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste("Time =", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot$add_text(
        "LPME",
        position = "upper_left",
        font_size = 12
      )
    }

    if (time_idx == length(time_points)) {
      adni_sil_plot$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    time_sil_plot$add_mesh(
      lpme_projection$boundary,
      show_edges = TRUE,
      color = "blue",
      line_width = line_width,
      opacity = opacity,
      label = "LPME"
    )

    adni_plot$subplot(2L, as.integer(time_idx - 1))
    adni_plot$add_mesh(
      pme_mesh$mesh,
      color = "green",
      show_edges = TRUE,
      opacity = 0.75,
      label = "PME"
    )

    if (time_idx == 1) {
      adni_plot$add_text(
        "PME",
        position = "lower_left",
        font_size = 12
      )
    }

    adni_sil_plot$subplot(1L, 0L)
    adni_sil_plot$add_mesh(
      pme_projection$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste("Time =", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot$add_text(
        "PME",
        position = "upper_left",
        font_size = 12
      )
    }

    if (time_idx == length(time_points)) {
      adni_sil_plot$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    time_sil_plot$add_mesh(
      pme_projection$boundary,
      show_edges = TRUE,
      color = "green",
      line_width = line_width,
      opacity = opacity,
      label = "PME"
    )

    adni_plot$subplot(3L, as.integer(time_idx - 1))
    adni_plot$add_mesh(
      pc_mesh$mesh,
      color = "purple",
      show_edges = TRUE,
      opacity = 0.75,
      label = "PS"
    )

    if (time_idx == 1) {
      adni_plot$add_text(
        "PS",
        position = "lower_left",
        font_size = 12
      )
    }

    adni_sil_plot$subplot(1L, 1L)
    adni_sil_plot$add_mesh(
      pc_projection$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste("Time =", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot$add_text(
        "PS",
        position = "upper_left",
        font_size = 12
      )
    }
    if (time_idx == length(time_points)) {
      adni_sil_plot$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    time_sil_plot$add_mesh(
      pc_projection$boundary,
      show_edges = TRUE,
      color = "purple",
      line_width = line_width,
      opacity = opacity,
      label = "PS"
    )

    time_sil_plot$add_legend(
      bcolor = "white",
      border = TRUE,
      loc = "upper right",
      face = "line"
    )
  }

  alg_labels <- dict(
    "Data" = "red",
    "LPME" = "blue",
    "PME" = "green",
    "PS" = "purple"
  )

  adni_plot$link_views()
  adni_plot$camera_position <- "yz"
  adni_plot$camera$azimuth <- -45
  # adni_plot$add_legend(
  #   alg_labels,
  #   bcolor = "white",
  #   border = TRUE,
  #   loc = "upper right",
  #   face = "line"
  # )

  adni_sil_plot$link_views()
  adni_sil_plot$camera_position <- "xy"
  # adni_sil_plot$add_scalar_bar(cmap = "Spectral", clim = c(0, 1))

  time_sil_plot$link_views()
  time_sil_plot$camera_position <- "xy"
  # adni_plot$add_legend(
  #   alg_labels,
  #   bcolor = "white",
  #   border = TRUE,
  #   loc = "upper right",
  #   face = "line"
  # )

  list(
    mesh_plot = adni_plot,
    sil_plot = adni_sil_plot,
    time_sil_plot = time_sil_plot
  )
}
