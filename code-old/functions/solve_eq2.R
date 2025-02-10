solve_eq2 <- function(E, T_val, Y, I, w, d, D) {
  M1 <- cbind(E + (w * diag(rep(1, I))), T_val)
  M2 <- cbind(t(T_val), matrix(0, ncol = d + 1, nrow = d + 1))
  M <- rbind(M1, M2)
  b <- rbind(Y, matrix(0, nrow = d + 1, ncol = D))
  sol <- ginv(M) %*% b
  return(sol)
}
