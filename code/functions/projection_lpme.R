projection_lpme <- function(x, f, initial.guess) {
  nlm_est <- try(
    nlm(
      function(t) dist_euclideanC(x = x, f(matrix(c(initial.guess[1], t), nrow = 1))),
      p = initial.guess[-1]
    )
  )
  if (class(nlm_est) == "try-error") {
    opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-07)
    nlopt_est <- try(
      nloptr(
        x0 = initial.guess[-1],
        function(t) dist_euclideanC(x = x, f(c(initial.guess[1], t))),
        opts = opts
      )
    )
    if (class(nlopt_est) == "try-error") {
      return(NULL)
    } else {
      return(c(initial.guess[1], nlopt_est$solution))
    }
  } else {
    return(c(initial.guess[1], nlm_est$estimate))
  }
}
