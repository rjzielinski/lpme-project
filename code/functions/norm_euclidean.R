norm_euclidean <- function(x) {
  # Norm function in a Euclidean space of any dimension
  norm_val <- sum(as.vector(x) ^ 2) %>%
    sqrt()
  return(norm_val)
}