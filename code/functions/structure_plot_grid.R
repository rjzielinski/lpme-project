structure_plot_grid <- function(
  adni_est,
  alpha_vals = exp(seq(-3, 0, 0.25)),
  opacity = 0.75,
  line_width = 2
) {
  require("reticulate")
  use_condaenv("lpme")
  pv <- import("pyvista")

  require(here)

  source(here("code/functions/estimate_mesh_volume_poisson.R"))
  source(here("code/functions/mesh_projection.R"))
  time_points <- unique(adni_est$data$data$time_from_bl)

  palette <- colorRampPalette(c(brewer.pal(8, "Spectral")))
  palette_colors <- palette(101)
  time_colors <- palette_colors[ceiling((time_points + 1e-10) * 100)]

  adni_plot <- pv$Plotter(shape = tuple(4L, as.integer(length(time_points))))
  adni_sil_plot1 <- pv$Plotter(shape = tuple(2L, 2L))
  adni_sil_plot2 <- pv$Plotter(shape = tuple(2L, 2L))
  adni_sil_plot3 <- pv$Plotter(shape = tuple(2L, 2L))

  time_sil_plot1 <- pv$Plotter(
    shape = tuple(2L, as.integer(ceiling(length(time_points) / 2)))
  )
  time_sil_plot2 <- pv$Plotter(
    shape = tuple(2L, as.integer(ceiling(length(time_points) / 2)))
  )
  time_sil_plot3 <- pv$Plotter(
    shape = tuple(2L, as.integer(ceiling(length(time_points) / 2)))
  )

  for (time_idx in seq_along(time_points)) {
    temp_data <- adni_est$data$data |>
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

    data_mesh <- estimate_mesh_volume_poisson(
      temp_data,
      alpha_vals = alpha_vals
    )
    lpme_mesh <- estimate_mesh_volume_poisson(
      temp_lpme_part_recon,
      alpha_vals = alpha_vals
    )
    pme_mesh <- estimate_mesh_volume_poisson(
      temp_pme_part_recon,
      alpha_vals = alpha_vals
    )
    pc_mesh <- estimate_mesh_volume_poisson(
      temp_pc_part_recon,
      alpha_vals = alpha_vals
    )

    data_projection1 <- mesh_projection(temp_data, axis = 1)
    data_projection2 <- mesh_projection(temp_data, axis = 2)
    data_projection3 <- mesh_projection(temp_data, axis = 3)

    lpme_projection1 <- mesh_projection(temp_lpme_part_recon, axis = 1)
    lpme_projection2 <- mesh_projection(temp_lpme_part_recon, axis = 2)
    lpme_projection3 <- mesh_projection(temp_lpme_part_recon, axis = 3)

    pme_projection1 <- mesh_projection(temp_pme_part_recon, axis = 1)
    pme_projection2 <- mesh_projection(temp_pme_part_recon, axis = 2)
    pme_projection3 <- mesh_projection(temp_pme_part_recon, axis = 3)

    pc_projection1 <- mesh_projection(temp_pc_part_recon, axis = 1)
    pc_projection2 <- mesh_projection(temp_pc_part_recon, axis = 2)
    pc_projection3 <- mesh_projection(temp_pc_part_recon, axis = 3)

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

    adni_sil_plot1$subplot(0L, 0L)
    adni_sil_plot1$add_mesh(
      data_projection1$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste0("Time = ", round(time_points[time_idx], 2))
    )

    adni_sil_plot2$subplot(0L, 0L)
    adni_sil_plot2$add_mesh(
      data_projection2$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste0("Time = ", round(time_points[time_idx], 2))
    )
    adni_sil_plot3$subplot(0L, 0L)
    adni_sil_plot3$add_mesh(
      data_projection3$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste0("Time = ", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot1$add_text(
        "Data",
        position = "upper_left",
        font_size = 12
      )
      adni_sil_plot2$add_text(
        "Data",
        position = "upper_left",
        font_size = 12
      )
      adni_sil_plot3$add_text(
        "Data",
        position = "upper_left",
        font_size = 12
      )
    }
    if (time_idx == length(time_points)) {
      adni_sil_plot1$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
      adni_sil_plot2$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
      adni_sil_plot3$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    time_sil_plot1$subplot(
      as.integer(round(time_idx / length(time_points), 0)),
      as.integer((time_idx - 1) %% ceiling(length(time_points) / 2))
    )
    time_sil_plot1$add_mesh(
      data_projection1$boundary,
      show_edges = TRUE,
      color = "red",
      line_width = line_width,
      opacity = opacity,
      label = "Data"
    )
    time_sil_plot1$add_text(
      paste("Time =", round(time_points[time_idx], 2)),
      position = "upper_left",
      font_size = 12
    )

    time_sil_plot2$subplot(
      as.integer(round(time_idx / length(time_points), 0)),
      as.integer((time_idx - 1) %% ceiling(length(time_points) / 2))
    )
    time_sil_plot2$add_mesh(
      data_projection2$boundary,
      show_edges = TRUE,
      color = "red",
      line_width = line_width,
      opacity = opacity,
      label = "Data"
    )
    time_sil_plot2$add_text(
      paste("Time =", round(time_points[time_idx], 2)),
      position = "upper_left",
      font_size = 12
    )

    time_sil_plot3$subplot(
      as.integer(round(time_idx / length(time_points), 0)),
      as.integer((time_idx - 1) %% ceiling(length(time_points) / 2))
    )
    time_sil_plot3$add_mesh(
      data_projection3$boundary,
      show_edges = TRUE,
      color = "red",
      line_width = line_width,
      opacity = opacity,
      label = "Data"
    )
    time_sil_plot3$add_text(
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

    adni_sil_plot1$subplot(0L, 1L)
    adni_sil_plot1$add_mesh(
      lpme_projection1$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste("Time =", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot1$add_text(
        "LPME",
        position = "upper_left",
        font_size = 12
      )
    }

    if (time_idx == length(time_points)) {
      adni_sil_plot1$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    adni_sil_plot2$subplot(0L, 1L)
    adni_sil_plot2$add_mesh(
      lpme_projection2$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste("Time =", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot2$add_text(
        "LPME",
        position = "upper_left",
        font_size = 12
      )
    }

    if (time_idx == length(time_points)) {
      adni_sil_plot2$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    adni_sil_plot3$subplot(0L, 1L)
    adni_sil_plot3$add_mesh(
      lpme_projection3$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste("Time =", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot3$add_text(
        "LPME",
        position = "upper_left",
        font_size = 12
      )
    }

    if (time_idx == length(time_points)) {
      adni_sil_plot3$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    time_sil_plot1$add_mesh(
      lpme_projection1$boundary,
      show_edges = TRUE,
      color = "blue",
      line_width = line_width,
      opacity = opacity,
      label = "LPME"
    )
    time_sil_plot2$add_mesh(
      lpme_projection2$boundary,
      show_edges = TRUE,
      color = "blue",
      line_width = line_width,
      opacity = opacity,
      label = "LPME"
    )
    time_sil_plot3$add_mesh(
      lpme_projection3$boundary,
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

    adni_sil_plot1$subplot(1L, 0L)
    adni_sil_plot1$add_mesh(
      pme_projection1$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste("Time =", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot1$add_text(
        "PME",
        position = "upper_left",
        font_size = 12
      )
    }

    if (time_idx == length(time_points)) {
      adni_sil_plot1$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    adni_sil_plot2$subplot(1L, 0L)
    adni_sil_plot2$add_mesh(
      pme_projection2$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste("Time =", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot2$add_text(
        "PME",
        position = "upper_left",
        font_size = 12
      )
    }

    if (time_idx == length(time_points)) {
      adni_sil_plot2$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    adni_sil_plot3$subplot(1L, 0L)
    adni_sil_plot3$add_mesh(
      pme_projection3$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste("Time =", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot3$add_text(
        "PME",
        position = "upper_left",
        font_size = 12
      )
    }

    if (time_idx == length(time_points)) {
      adni_sil_plot3$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    time_sil_plot1$add_mesh(
      pme_projection1$boundary,
      show_edges = TRUE,
      color = "green",
      line_width = line_width,
      opacity = opacity,
      label = "PME"
    )

    time_sil_plot2$add_mesh(
      pme_projection2$boundary,
      show_edges = TRUE,
      color = "green",
      line_width = line_width,
      opacity = opacity,
      label = "PME"
    )

    time_sil_plot3$add_mesh(
      pme_projection3$boundary,
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

    adni_sil_plot1$subplot(1L, 1L)
    adni_sil_plot1$add_mesh(
      pc_projection1$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste("Time =", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot1$add_text(
        "PS",
        position = "upper_left",
        font_size = 12
      )
    }
    if (time_idx == length(time_points)) {
      adni_sil_plot1$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    adni_sil_plot2$subplot(1L, 1L)
    adni_sil_plot2$add_mesh(
      pc_projection2$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste("Time =", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot2$add_text(
        "PS",
        position = "upper_left",
        font_size = 12
      )
    }
    if (time_idx == length(time_points)) {
      adni_sil_plot2$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    adni_sil_plot3$subplot(1L, 1L)
    adni_sil_plot3$add_mesh(
      pc_projection3$boundary,
      show_edges = TRUE,
      color = time_colors[time_idx],
      line_width = line_width,
      opacity = opacity,
      label = paste("Time =", round(time_points[time_idx], 2))
    )

    if (time_idx == 1) {
      adni_sil_plot3$add_text(
        "PS",
        position = "upper_left",
        font_size = 12
      )
    }
    if (time_idx == length(time_points)) {
      adni_sil_plot3$add_legend(
        bcolor = "white",
        border = TRUE,
        loc = "upper right",
        face = "line"
      )
    }

    time_sil_plot1$add_mesh(
      pc_projection1$boundary,
      show_edges = TRUE,
      color = "purple",
      line_width = line_width,
      opacity = opacity,
      label = "PS"
    )

    time_sil_plot1$add_legend(
      bcolor = "white",
      border = TRUE,
      loc = "upper right",
      face = "line"
    )

    time_sil_plot2$add_mesh(
      pc_projection2$boundary,
      show_edges = TRUE,
      color = "purple",
      line_width = line_width,
      opacity = opacity,
      label = "PS"
    )

    time_sil_plot2$add_legend(
      bcolor = "white",
      border = TRUE,
      loc = "upper right",
      face = "line"
    )

    time_sil_plot3$add_mesh(
      pc_projection3$boundary,
      show_edges = TRUE,
      color = "purple",
      line_width = line_width,
      opacity = opacity,
      label = "PS"
    )

    time_sil_plot3$add_legend(
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

  adni_sil_plot1$link_views()
  adni_sil_plot1$camera_position <- "yz"

  adni_sil_plot2$link_views()
  adni_sil_plot2$camera_position <- "xz"

  adni_sil_plot3$link_views()
  adni_sil_plot3$camera_position <- "xy"

  time_sil_plot1$link_views()
  time_sil_plot1$camera_position <- "yz"

  time_sil_plot2$link_views()
  time_sil_plot2$camera_position <- "xz"

  time_sil_plot3$link_views()
  time_sil_plot3$camera_position <- "xy"

  list(
    mesh_plot = adni_plot,
    sil_plots = list(
      dim1 = adni_sil_plot1,
      dim2 = adni_sil_plot2,
      dim3 = adni_sil_plot3
    ),
    time_sil_plots = list(
      dim1 = time_sil_plot1,
      dim2 = time_sil_plot2,
      dim3 = time_sil_plot3
    )
  )
}
