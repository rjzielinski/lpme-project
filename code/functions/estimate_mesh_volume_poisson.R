estimate_mesh_volume_poisson <- function(
  points,
  depth_val = 3,
  alpha_vals = seq(0.2, 2.0, by = 0.1),
  ball_radii = c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1),
  threshold = 0.005,
  plot_mesh = FALSE
) {
  require(dplyr, warn.conflicts = FALSE, quietly = TRUE)

  pcd <- o3d$geometry$PointCloud()
  pcd$points <- o3d$utility$Vector3dVector(points)

  pcd$estimate_normals()
  pcd$orient_normals_consistent_tangent_plane(100L)

  poisson_results <- o3d$geometry$TriangleMesh$create_from_point_cloud_poisson(
    pcd,
    depth = as.integer(depth_val),
    width = 0,
    scale = 1.1,
    linear_fit = TRUE,
    n_threads = 1L
  )

  ball_mesh <- o3d$geometry$TriangleMesh$create_from_point_cloud_ball_pivoting(
    pcd,
    o3d$utility$DoubleVector(ball_radii)
  )

  poisson_mesh <- poisson_results[[1]]
  poisson_densities <- np$array(poisson_results[[2]])

  if (poisson_mesh$is_watertight() == TRUE) {
    vertices <- np$asarray(poisson_mesh$vertices)
    faces <- np$asarray(poisson_mesh$triangles)

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
    pv_mesh <- pv$PolyData(vertices, faces_aug)
  } else if (ball_mesh$is_watertight() == TRUE) {
    vertices <- np$asarray(ball_mesh$vertices)
    faces <- np$asarray(ball_mesh$triangles)

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
    pv_mesh <- pv$PolyData(vertices, faces_aug)
  } else {
    # use Delaunay triangulation to estimate volume of point cloud
    # loop through possible alpha values to find stable volume estimate
    # assume that we have already loaded pyvista as pv through reticulate package

    cloud <- pv$PolyData(points)

    pv_mesh <- cloud$reconstruct_surface()

    #   mesh_list <- list()
    #   mesh_volumes <- vector(mode = "numeric", length = length(alpha_vals))
    #
    #   for (alpha_idx in seq_along(alpha_vals)) {
    #     mesh_list[[alpha_idx]] <- cloud$delaunay_3d(alpha = alpha_vals[alpha_idx])
    #     mesh_volumes[alpha_idx] <- mesh_list[[alpha_idx]]$volume
    #   }
    #
    #   mesh_vol_change <- c(
    #     NA,
    #     ((lead(mesh_volumes) - mesh_volumes) / mesh_volumes)[
    #       -length(mesh_volumes)
    #     ]
    #   )
    #
    #   n_elig <- length(which(abs(mesh_vol_change) < threshold))
    #
    #   first_idx <- ifelse(
    #     length(which(abs(mesh_vol_change) < threshold)) == 0,
    #     which.min(abs(mesh_vol_change)),
    #     min(which(abs(mesh_vol_change) < threshold))
    #   )
    #
    #   out_mesh <- mesh_list[[first_idx]]
  }

  # mf <- pymeshfix$MeshFix(pv_mesh)
  # mf$repair(verbose = TRUE)
  # out_mesh <- mf$mesh
  #
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
