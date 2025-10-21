pv_to_o3d <- function(pv_mesh) {
  require(reticulate, warn.conflicts = FALSE, quietly = TRUE)
  use_condaenv("lpme")
  np <- import("numpy")
  o3d <- import("open3d")
  mf <- import("pymeshfix")
  plt <- import("matplotlib.pyplot")
  pv <- import("pyvista")

  pv_mesh <- pv_mesh$triangulate()
  vertices <- np$asarray(pv_mesh$points)
  faces <- np$asarray(pv_mesh$faces)
  faces <- array_reshape(faces, c(length(faces) / 4, 4))[, -1]

  o3d_mesh <- o3d$geometry$TriangleMesh()
  o3d_mesh$vertices <- o3d$utility$Vector3dVector(vertices)
  o3d_mesh$triangles <- o3d$utility$Vector3iVector(faces)

  o3d_mesh$compute_vertex_normals()
  o3d_mesh$compute_triangle_normals()

  return(o3d_mesh)
}

o3d_to_pv <- function(o3d_mesh) {
  require(reticulate, warn.conflicts = FALSE, quietly = TRUE)
  use_condaenv("lpme")
  np <- import("numpy")
  o3d <- import("open3d")
  mf <- import("pymeshfix")
  plt <- import("matplotlib.pyplot")
  pv <- import("pyvista")

  vertices <- np$asarray(o3d_mesh$vertices)
  faces <- np$asarray(o3d_mesh$triangles)

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

  return(pv_mesh)
}
