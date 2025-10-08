estimate_mesh_volume_poisson <- function(
  points,
  depth_val = 3,
  alpha_vals = seq(0.2, 2.0, by = 0.1),
  threshold = 0.005,
  plot_mesh = FALSE
) {
  require(dplyr, warn.conflicts = FALSE, quietly = TRUE)

  pcd <- o3d$geometry$PointCloud()
  pcd$points <- o3d$utility$Vector3dVector(points)

  pcd$estimate_normals()
  pcd$orient_normals_consistent_tangent_plane(100L)

  results <- o3d$geometry$TriangleMesh$create_from_point_cloud_poisson(
    pcd,
    depth = as.integer(depth_val),
    width = 0,
    scale = 1.1,
    linear_fit = TRUE,
    n_threads = 1L
  )

  mesh <- results[[1]]
  densities <- np$array(results[[2]])

  if (mesh$is_watertight() == TRUE) {
    vertices <- np$asarray(mesh$vertices)
    faces <- np$asarray(mesh$triangles)

    # PyVista requires faces in the format [3, p1, p2, p3, 3, p4, p5, p6, ...]
    # We'll use numpy's functions to prepend a '3' to each face

    # Get the number of faces (triangles)
    num_faces <- dim(faces)[1]

    # Create a 1D NumPy array of '3's, one for each face
    int_vec <- np$full(as.integer(num_faces), 3L)

    # Combine the '3's column with the face indices column-wise
    # np$column_stack is the function equivalent of Python's np.c_
    faces_aug <- cbind(int_vec, faces)
    mode(faces_aug) <- "integer"

    # Create and return the PyVista PolyData object
    out_mesh <- pv$PolyData(vertices, faces_aug)
  } else {
    # use Delaunay triangulation to estimate volume of point cloud
    # loop through possible alpha values to find stable volume estimate
    # assume that we have already loaded pyvista as pv through reticulate package

    cloud <- pv$PolyData(points)

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
    out_mesh <- mesh_list[[first_idx]]
  }

  # densities_r <- np$array(densities)
  # density_threshold <- quantile(densities_r, threshold)

  # vertices_to_remove <- densities_r < density_threshold
  # mesh$remove_vertices_by_mask(vertices_to_remove)

  if (plot_mesh == TRUE) {
    mesh_plot <- pv$Plotter()
    mesh_plot$add_mesh(out_mesh)
    mesh_plot$show()
  }

  out_list <- list(
    mesh = out_mesh,
    volume = out_mesh$volume
  )

  return(out_list)
}
