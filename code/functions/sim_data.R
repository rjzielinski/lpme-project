sim_data <- function(time_val, case, noise, amp_noise, period_noise, N = 1000) {
  manifolds <- list(
    function(tau, amp_noise, period_noise) {
      return(c(tau, amp_noise[1] * sin(period_noise[1] * tau + pi / 2)))
    },
    function(tau, amp_noise, period_noise) {
      return(c(tau, amp_noise[1] * sin(period_noise[1] * tau)))
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          amp_noise[1] * cos(period_noise[1] * tau),
          amp_noise[2] * sin(period_noise[2] * tau)
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          amp_noise[1] * cos(period_noise[1] * tau),
          amp_noise[2] * sin(period_noise[2] * tau)
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          tau,
          (tau * amp_noise[1] + period_noise[1]) ^ 2,
          (tau * amp_noise[2] + period_noise[2]) ^ 3
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          tau,
          amp_noise[1] * cos(period_noise[1] * tau),
          amp_noise[2] * sin(period_noise[2] * tau)
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          period_noise[1] * tau[1],
          period_noise[2] * tau[2],
          amp_noise[1] * (amp_noise[2] * norm_euclidean(period_noise[1:d] * tau) ^ 2)
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          (period_noise[1] * amp_noise[1] * tau[1]) * cos(amp_noise[1] * tau[1]),
          (period_noise[2] * amp_noise[2] * tau[1]) * sin(amp_noise[2] * tau[1]),
          tau[2]
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      return(
        c(
          amp_noise[1] * sin(period_noise[1] * tau[1]) * cos(period_noise[2] * tau[2]),
          amp_noise[1] * sin(period_noise[1] * tau[1]) * sin(period_noise[2] * tau[2]),
          amp_noise[1] * cos(period_noise[1] * tau[1])
        )
      )
    },
    function(tau, amp_noise, period_noise) {
      r <- 1 + amp_noise[1] * (tau[1] + 1) * sqrt(tau[2] + 1)
      return(
        c(
          r * sin(period_noise[1] * tau[1]) * cos(period_noise[2] * tau[2]),
          r * sin(period_noise[1] * tau[1]) * sin(period_noise[2] * tau[2]),
          r * cos(period_noise[1] * tau[1])
        )
      )
    }
  )

  D <- case_when(
    case == 1 ~ 2,
    case == 2 ~ 2,
    case == 3 ~ 2,
    case == 4 ~ 2,
    case == 5 ~ 3,
    case == 6 ~ 3,
    case == 7 ~ 3,
    case == 8 ~ 3,
    case == 9 ~ 3,
    case == 10 ~ 3
  )
  d <- case_when(
    case == 1 ~ 1,
    case == 2 ~ 1,
    case == 3 ~ 1,
    case == 4 ~ 1,
    case == 5 ~ 1,
    case == 6 ~ 1,
    case == 7 ~ 2,
    case == 8 ~ 2,
    case == 9 ~ 2,
    case == 10 ~ 2
  )

  manifold <- manifolds[[case]]

  I <- N
  t <- matrix(NA, nrow = I, ncol = d)
  X <- matrix(NA, nrow = I, ncol = D)
  noise_vals <- rnorm(I * D, mean = 0, sd = noise) %>%
    matrix(nrow = I, ncol = D)

  if (case == 1) {
    t[, 1] <- rnorm(I, mean = 0, sd = 1)
  } else if (case == 2) {
    t[, 1] <- runif(I, min = -3 * pi, max = 3 * pi)
  } else if (case == 3) {
    t[, 1] <- runif(I, min = -0.8 * pi, max = 0.5 * pi)
  } else if (case == 4) {
    t[, 1] <- runif(I, time_val * (pi / 4), (time_val + 5) * (pi / 4))
  } else if (case == 5) {
    t[, 1] <- runif(I, min = -1, max = 1)
  } else if (case == 6) {
    t[, 1] <- runif(I, min = 0, max = 3 * pi)
  } else if (case == 7) {
    t[, 1] <- runif(I, min = -1, max = 1)
    t[, 2] <- runif(I, min = -1, max = 1)
  } else if (case == 8) {
    t[, 1] <- runif(I, min = 0, max = 10)
    t[, 2] <- runif(I, min = -1, max = 1)
  } else if (case == 9) {
    t[, 1] <- runif(I, min = 0, max = pi)
    t[, 2] <- runif(I, min = 0, max = 2 * pi)
  } else if (case == 10) {
    t[, 1] <- runif(I, min = 0, max = pi)
    t[, 2] <- runif(I, min = 0, max = 2 * pi)
  }

  amp_mean <- 1
  per_mean <- 1

  # if (case %in% 1:10) {
  #   amp_mean <- 1
  #   per_mean <- 1
  # }

  amp_noise <- rnorm(D, mean = amp_mean, sd = amp_noise)
  period_noise <- rnorm(D, mean = per_mean, sd = period_noise)

  X <- map(
    1:nrow(t),
    ~ manifold(
      t[.x, ],
      amp_noise,
      period_noise
    )
  ) %>%
    unlist() %>%
    matrix(ncol = D, byrow = TRUE)

  data.points <- X + noise_vals
  data.points <- cbind(time_val, data.points)
  return(data.points)
}
