estimate_mesh_volume_poisson <- function(
  points,
  depth_val = 3,
  alpha_vals = exp(seq(-3, 0, 0.25)),
  ball_radii = c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1),
  threshold = 0.01,
  mesh_error_vol = 0.5,
  fix_meshes = FALSE,
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

  # Create and return the PyVista PolyData object
  init_poisson_vol <- ifelse(
    poisson_mesh$is_watertight(),
    poisson_mesh$get_volume(),
    0
  )
  pv_mesh <- o3d_to_pv(poisson_mesh)
  mesh_volume <- init_poisson_vol

  meshfix <- mf$MeshFix(pv_mesh)
  holes <- meshfix$extract_holes()

  if (poisson_mesh$is_watertight() == FALSE && fix_meshes == TRUE) {
    meshfix$repair()
    pv_mesh <- meshfix$mesh
    poisson_mesh <- pv_to_o3d(pv_mesh)
  }

  if (
    poisson_mesh$is_watertight() == FALSE || init_poisson_vol < mesh_error_vol
  ) {
    # if there are still holes, switch to ball-pivoting mesh

    ball_mesh <- o3d$geometry$TriangleMesh$create_from_point_cloud_ball_pivoting(
      pcd,
      o3d$utility$DoubleVector(ball_radii)
    )

    pv_mesh <- o3d_to_pv(ball_mesh)
    init_ball_vol <- ifelse(
      ball_mesh$is_watertight(),
      ball_mesh$get_volume(),
      0
    )
    mesh_volume <- init_ball_vol

    meshfix <- mf$MeshFix(pv_mesh)
    holes <- meshfix$extract_holes()

    if (ball_mesh$is_watertight() == FALSE && fix_meshes == TRUE) {
      meshfix$repair()

      pv_mesh <- meshfix$mesh
      ball_mesh <- pv_to_o3d(pv_mesh)
    }

    if (ball_mesh$is_watertight() == FALSE || init_ball_vol < mesh_error_vol) {
      # use alpha shapes to estimate volume of point cloud
      # loop through possible alpha values to find stable volume estimate

      mesh_list <- list()
      mesh_volumes <- vector(mode = "numeric", length = length(alpha_vals))

      watertight <- FALSE
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
          mesh_volumes[alpha_idx] <- ifelse(
            alpha_mesh$is_watertight(),
            alpha_mesh$get_volume(),
            NA
          )

          meshfix <- mf$MeshFix(pv_mesh)
          holes <- meshfix$extract_holes()

          if (alpha_mesh$is_watertight() == TRUE) {
            if (alpha_mesh$get_volume() > mesh_error_vol) {
              pv_mesh <- o3d_to_pv(alpha_mesh)
              mesh_volume <- alpha_mesh$get_volume()

              watertight <- TRUE
              break
            }
          }
        } else {
          mesh_list[[alpha_idx]] <- NULL
          mesh_volumes[alpha_idx] <- NA
        }
      }

      if (watertight == FALSE) {
        convex_hull_trimesh <- pcd$compute_convex_hull()[[1]]
        pv_mesh <- o3d_to_pv(convex_hull_trimesh)
        mesh_volume <- convex_hull_trimesh$get_volume()
      }
    }
  }

  out_mesh <- pv_mesh

  if (plot_mesh == TRUE) {
    mesh_plot <- pv$Plotter()
    mesh_plot$add_mesh(out_mesh)
    mesh_plot$show()
  }

  out_list <- list(
    mesh = out_mesh,
    volume = mesh_volume
  )

  return(out_list)
}
