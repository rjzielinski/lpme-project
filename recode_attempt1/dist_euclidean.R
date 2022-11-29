dist.euclidean <- function(x, y) {
  # Distance function in a Euclidean space of any dimension
  return(norm.euclidean(x - y))
}