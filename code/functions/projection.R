projection <- function(x, f, initial.guess) {
  est <- nlm(function(t) dist_euclideanC(x = x, f(t)), p = initial.guess)
  # opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-04)
  # est <- nloptr(
  #   x0 = initial.guess,
  #   function(t) dist_euclideanC(x = x, f(t)),
  #   opts = opts
  # )
  return(est$estimate)
}
