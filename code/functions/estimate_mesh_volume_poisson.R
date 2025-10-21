estimate_mesh_volume_poisson <- function(
  points,
  depth_val = 3,
  alpha_vals = exp(seq(-2, 0, 0.25)),
  ball_radii = c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1),
  threshold = 0.01,
  mesh_error_vol = 0.5,
  plot_mesh = FALSE
) {
  require(dplyr, warn.conflicts = FALSE, quietly = TRUE)
  require(here)

  require(reticulate, warn.conflicts = FALSE, quietly = TRUE)
  use_condaenv("lpme")
  np <- import("numpy")
  o3d <- import("open3d")
  mf <- import("pymeshfix")
  plt <- import("matplotlib.pyplot")
  pv <- import("pyvista")
  warnings <- import("warnings")

  warnings$simplefilter("ignore")

  source(here("code/functions/convert_mesh_lib.R"))

  pcd <- o3d$geometry$PointCloud()
  pcd$points <- o3d$utility$Vector3dVector(points)

  pcd$estimate_normals()
  pcd$orient_normals_consistent_tangent_plane(20L)

  poisson_results <- o3d$geometry$TriangleMesh$create_from_point_cloud_poisson(
    pcd,
    depth = as.integer(depth_val),
    width = 0,
    scale = 1.1,
    linear_fit = TRUE,
    n_threads = 1L
  )
  poisson_mesh <- poisson_results[[1]]
  poisson_densities <- np$array(poisson_results[[2]])

  poisson_mesh <- poisson_mesh$remove_degenerate_triangles()
  poisson_mesh <- poisson_mesh$remove_duplicated_triangles()
  poisson_mesh <- poisson_mesh$remove_non_manifold_edges()
  poisson_mesh <- poisson_mesh$remove_unreferenced_vertices()

  poisson_mesh_filled <- o3d$t$geometry$TriangleMesh$from_legacy(poisson_mesh)
  poisson_mesh_filled <- poisson_mesh_filled$fill_holes(1000000)

  poisson_mesh <- poisson_mesh_filled$to_legacy()

  # Create and return the PyVista PolyData object
  pv_mesh <- o3d_to_pv(poisson_mesh)
  init_poisson_vol <- pv_mesh$volume

  meshfix <- mf$MeshFix(pv_mesh)
  holes <- meshfix$extract_holes()

  if (holes$n_cells > 0 || init_poisson_vol < mesh_error_vol) {
    # if there are still holes, switch to ball-pivoting mesh
    ball_mesh <- o3d$geometry$TriangleMesh$create_from_point_cloud_ball_pivoting(
      pcd,
      o3d$utility$DoubleVector(ball_radii)
    )

    ball_mesh <- ball_mesh$remove_degenerate_triangles()
    ball_mesh <- ball_mesh$remove_duplicated_triangles()
    ball_mesh <- ball_mesh$remove_non_manifold_edges()
    ball_mesh <- ball_mesh$remove_unreferenced_vertices()

    ball_mesh_filled <- o3d$t$geometry$TriangleMesh$from_legacy(ball_mesh)
    ball_mesh_filled <- ball_mesh_filled$fill_holes(1000000)

    ball_mesh <- ball_mesh_filled$to_legacy()

    pv_mesh <- o3d_to_pv(ball_mesh)
    init_ball_vol <- pv_mesh$volume

    # check ball mesh for holes
    meshfix <- mf$MeshFix(pv_mesh)
    ball_holes <- meshfix$extract_holes()

    if (ball_holes$n_cells > 0 || init_ball_vol < mesh_error_vol) {
      # use alpha shapes to estimate volume of point cloud
      # loop through possible alpha values to find stable volume estimate

      mesh_list <- list()
      mesh_volumes <- vector(mode = "numeric", length = length(alpha_vals))

      convex_hull <- o3d$geometry$TetraMesh$create_from_point_cloud(pcd)
      for (alpha_idx in seq_along(alpha_vals)) {
        alpha_mesh <- o3d$geometry$TriangleMesh$create_from_point_cloud_alpha_shape(
          pcd,
          alpha_vals[alpha_idx],
          convex_hull[[1]],
          convex_hull[[2]]
        )

        # alpha_mesh$compute_vertex_normals()

        if (!alpha_mesh$is_empty()) {
          mesh_list[[alpha_idx]] <- o3d_to_pv(alpha_mesh)
          mesh_volumes[alpha_idx] <- mesh_list[[alpha_idx]]$volume
        } else {
          mesh_list[[alpha_idx]] <- NULL
          mesh_volumes[alpha_idx] <- NA
        }

        if (alpha_mesh$is_watertight() == TRUE) {
          if (alpha_mesh$get_volume() > mesh_error_vol) {
            pv_mesh <- o3d_to_pv(alpha_mesh)
            break
          }
        }

        # o3d$visualization$draw_geometries(
        #   list(alpha_mesh),
        #   mesh_show_back_face = TRUE
        # )

        # mesh_list[[alpha_idx]] <- o3d_to_pv(alpha_mesh)
        # mesh_volumes[alpha_idx] <- mesh_list[[alpha_idx]]$volume
        #
        # mesh_vol_change <- c(
        #   NA,
        #   ((lead(mesh_volumes) - mesh_volumes) / mesh_volumes)[
        #     -length(mesh_volumes)
        #   ]
        # )
      }

      # cloud <- pv$PolyData(points)
      #
      # mesh_list <- list()
      # mesh_volumes <- vector(mode = "numeric", length = length(alpha_vals))
      #
      # for (alpha_idx in seq_along(alpha_vals)) {
      #   alpha_mesh <- cloud$delaunay_3d(
      #     alpha = alpha_vals[alpha_idx]
      #   )$extract_surface()
      #
      #   mesh_list[[alpha_idx]] <- alpha_mesh
      #   mesh_volumes[alpha_idx] <- mesh_list[[alpha_idx]]$volume
      #
      #   mesh_vol_change <- c(
      #     NA,
      #     ((lead(mesh_volumes) - mesh_volumes) / mesh_volumes)[
      #       -length(mesh_volumes)
      #     ]
      #   )
      #
      #   if (sum(abs(mesh_vol_change) < threshold, na.rm = TRUE) > 0) {
      #     pv_mesh <- mesh_list[[alpha_idx]]
      #     break
      #   }
      # }

      mesh_vol_change <- c(
        NA,
        ((lead(mesh_volumes) - mesh_volumes) / mesh_volumes)[
          -length(mesh_volumes)
        ]
      )

      n_elig <- length(which(abs(mesh_vol_change) < threshold))

      if (n_elig == 0) {
        pv_mesh <- mesh_list[[which.min(abs(mesh_vol_change))]]
      }
    }
  }

  out_mesh <- pv_mesh
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
