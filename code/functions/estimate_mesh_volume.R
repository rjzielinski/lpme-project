estimate_mesh_volume <- function(
  cloud,
  alpha_vals,
  threshold = 0.001,
  plot_mesh = FALSE
) {
  points <- temp_data
  pcd <- o3d$geometry$PointCloud()
  pcd$points <- o3d$utility$Vector3dVector(points)

  pcd$estimate_normals()
  pcd$orient_normals_consistent_tangent_plane(100L)

  depth <- 3

  results <- o3d$geometry$TriangleMesh$create_from_point_cloud_poisson(
    pcd,
    depth = as.integer(depth),
    width = 0,
    scale = 1.1,
    linear_fit = FALSE
  )

  mesh <- results[[1]]
  densities <- results[[2]]

  densities_r <- np$array(densities)
  # density_threshold <- quantile(densities_r, 0.025)

  vertices_to_remove <- densities_r < density_threshold
  mesh$remove_vertices_by_mask(vertices_to_remove)

  mesh$paint_uniform_color(c(0.6, 0.8, 1.0))
  mesh$compute_vertex_normals()

  o3d$visualization$draw_geometries(list(mesh))

  # use Delaunay triangulation to estimate volume of point cloud
  # loop through possible alpha values to find stable volume estimate
  # assume that we have already loaded pyvista as pv through reticulate package

  require(dplyr, warn.conflicts = FALSE, quietly = TRUE)

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

  first_idx <- min(which(abs(mesh_vol_change) < threshold))
  # first_idx <- which.min(abs(mesh_vol_change))
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
