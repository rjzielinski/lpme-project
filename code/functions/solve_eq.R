solve_eq <- function(E, W, t_val, X, w, d, D) {
  M1 <- cbind(
    2 * E %*% W %*% E + 2 * w * E,
    2 * E %*% W %*% t_val,
    t_val
  )
  M2 <- cbind(
    2 * t(t_val) %*% W %*% E,
    2 * t(t_val) %*% W %*% t_val,
    matrix(0, ncol = d + 1, nrow = d + 1)
  )
  M3 <- cbind(
    t(t_val),
    matrix(0, ncol = d + 1, nrow = d + 1),
    matrix(0, ncol = d + 1, nrow = d + 1)
  )
  M <- rbind(M1, M2, M3) # The coefficient matrix of the linear equations

  # The nonhomogeneous term of the linear equations
  b <- rbind(
    2 * E %*% W %*% X,
    2 * t(t_val) %*% W %*% X,
    matrix(0, nrow = d + 1, ncol = D)
  )
  sol <- ginv(M) %*% b # Solve the linear equations
  return(sol)
}
