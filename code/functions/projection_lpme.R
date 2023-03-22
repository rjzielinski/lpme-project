projection_lpme <- function(x, f, initial.guess) {
  est <- nlm(function(t) dist_euclideanC(x = x, f(matrix(c(initial.guess[1], t), nrow = 1))), p = initial.guess[-1])
  # opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-04)
  # est <- nloptr(
  #   x0 = initial.guess,
  #   function(t) dist_euclideanC(x = x, f(t)),
  #   opts = opts
  # )
  return(c(initial.guess[1], est$estimate))
}
