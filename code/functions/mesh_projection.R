mesh_projection <- function(data) {
  require("reticulate")
  use_condaenv("lpme")
  pv <- import("pyvista")

  cloud <- pv$PolyData(data)

  origin = cloud$center
  origin[[length(origin)]] <- cloud$length / 3.0
  projected <- cloud$project_points_to_plane(origin = origin)

  out_mesh <- projected$delaunay_2d(alpha = 0.25)
  edges <- out_mesh$extract_feature_edges(
    boundary_edges = TRUE,
    feature_edges = FALSE,
    manifold_edges = FALSE,
    non_manifold_edges = FALSE
  )
  out_list <- list(
    mesh = out_mesh,
    boundary = edges
  )
}
