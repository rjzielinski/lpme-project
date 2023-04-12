projection_lpme <- function(x, f, initial.guess) {
  nlm_est <- try(
    nlm(
      function(t) dist_euclideanC(x = x, f(matrix(c(initial.guess[1], t), nrow = 1))),
      p = initial.guess[-1]
    )
  )
  count <- 1
  if (class(nlm_est) == "try-error") {
    while (class(nlm_est) == "try-error" & count < 10) {
      nlm_est <- try(
        nlm(
          function(t) dist_euclideanC(x = x, f(matrix(c(initial.guess[1], t), nrow = 1))),
          p = initial.guess[-1]
        )
      )
      count <- count + 1
    }
    if (class(nlm_est) != "try-error") {
      return(c(initial.guess[1], nlm_est$estimate))
    } else {
      opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-04)
      nlopt_est <- try(
        nloptr(
          x0 = initial.guess,
          function(t) dist_euclideanC(x = x, f(t)),
          opts = opts
        )
      )
      if (class(nlopt_est) == "try-error") {
        return(NA)
      } else {
        return(c(initial.guess[1], nlopt_est$solution))
      }
    }
  } else {
    return(c(initial.guess[1], nlm_est$estimate))
  }
}
