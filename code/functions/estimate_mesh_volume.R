estimate_mesh_volume <- function(
  cloud,
  alpha_vals,
  threshold = 0.005,
  plot_mesh = FALSE
) {
  # use Delaunay triangulation to estimate volume of point cloud
  # loop through possible alpha values to find stable volume estimate
  # assume that we have already loaded pyvista as pv through reticulate package

  mesh_list <- list()
  mesh_volumes <- vector(mode = "numeric", length = length(alpha_vals))

  for (alpha_idx in seq_along(alpha_vals)) {
    mesh_list[[alpha_idx]] <- cloud$delaunay_3d(alpha = alpha_vals[alpha_idx])
    mesh_volumes[alpha_idx] <- mesh_list[[alpha_idx]]$volume
  }

  mesh_vol_change <- c(
    NA,
    ((lead(mesh_volumes) - mesh_volumes) / mesh_volumes)[
      -length(mesh_volumes)
    ]
  )

  first_idx <- min(which(mesh_vol_change < threshold))
  out_mesh <- mesh_list[[first_idx]]
  out_vol <- out_mesh$volume

  if (plot_mesh == TRUE) {
    mesh_plot <- pv$Plotter()
    mesh_plot$add_mesh(out_mesh)
    mesh_plot$show()
  }

  out_list <- list(
    mesh = out_mesh,
    volume = out_vol,
    alpha = alpha_vals[first_idx]
  )
  return(out_list)
}
