norm.euclidean <- function(x) {
  # Norm function in a Euclidean space of any dimension
  require(Matrix)
  return(norm(matrix(x, ncol = 1), type = "F"))
}