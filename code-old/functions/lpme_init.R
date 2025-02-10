lpme_init <- function(df, d, alpha = 0.05, max.comp = 100, epsilon = 0.05, max.iter = 100, init = "first") {

  # init argument can take three values:
  # "first" indicates that an isomap initialization run on the first time point is used
  # "full" indicates that an isomap initialization run on the full dataset (all time points) is used
  # "separate" indicates that each time point will use its own isomap initialization

  source("code/pme.R")
  require(plotly)

  time_points <- df[, 1] %>%
    unique()

  pme_results <- list()
  funcs <- list()
  x_test <- list()
  r <- list()

  msd <- list()

  if (init %in% c("first", "full")) {
    if (init == "first") {
      init_df <- df[df[, 1] == time_points[1], -1]
    } else if (init == "full") {
      init_df <- df[, -1]
    }

    init_dimension.size <- dim(init_df)
    init_D <- init_dimension.size[2]
    init_n <- init_dimension.size[1]
    lambda <- 4 - d

    init_N0 <- 20 * init_D
    init_est <- hdmde(init_df, init_N0, alpha, max.comp)
    init_order <- order(init_est$mu[, 1])
    init_theta.hat <- init_est$theta.hat[init_order]
    init_centers <- init_est$mu[init_order, ]
    init_sigma <- init_est$sigma
    init_W <- diag(init_theta.hat)
    init_X <- init_est$mu[init_order, ]
    init_I <- length(init_theta.hat)

    init_dissimilarity.matrix <- as.matrix(dist(init_X))
    init_isomap <- isomap(init_dissimilarity.matrix, ndim = d, k = 10)
  }

  for (idx in 1:length(time_points)) {
    df_temp <- df[df[, 1] == time_points[idx], ]
    if (init %in% c("first", "full")) {
      pme_results[[idx]] <- pme(
        x.obs = df_temp[, -1],
        d = d,
        initialization = list(init_est, init_isomap)
      )
    } else {
      pme_results[[idx]] <- pme(
        x.obs = df_temp[, -1],
        d = d
      )
    }

    msd[[idx]] <- min(pme_results[[idx]]$MSD)
    funcs[[idx]] <- pme_results[[idx]]$embedding.map
    r_test <- seq(
      from = -10,
      to = 10,
      length.out = dim(pme_results[[idx]]$knots)[1]
    )
    r_list <- lapply(numeric(d), function(x) r_test)
    r_mat <- as.matrix(expand.grid(r_list))
    r_length <- dim(r_mat)[1]
    x_temp <- apply(
      r_mat,
      1,
      funcs[[idx]]
    ) %>%
      matrix(nrow = dim(r_mat)[1], byrow = TRUE)

    idx_inrange <- matrix(nrow = dim(x_temp)[1], ncol = dim(x_temp)[2])
    for (dim_idx in 1:dim(x_temp)[2]) {
      idx_range <- max(df_temp[, dim_idx + 1]) - min(df_temp[, dim_idx + 1])
      idx_min <- min(df_temp[, dim_idx + 1]) - (0.2 * idx_range)
      idx_max <- max(df_temp[, dim_idx + 1]) + (0.2 * idx_range)
      idx_inrange[, dim_idx] <- (x_temp[, dim_idx] > idx_min) & (x_temp[, dim_idx] < idx_max)
    }

    r_inrange <- rowSums(idx_inrange) == dim(x_temp)[2]
    if (sum(r_inrange) == 0) {
      r_inrange <- rowSums(idx_inrange) > 0
    }
    r_min <- min(r_mat[r_inrange, ])
    r_min <- ifelse(is.na(r_min), -10, r_min)
    r_max <- max(r_mat[r_inrange, ])
    r_max <- ifelse(is.na(r_max), 10, r_max)
    r_test <- seq(
      r_min,
      r_max,
      length.out = max(10, round((dim(pme_results[[idx]]$knots)[1]) ^ (1 / d)))
    )
    r_list <- lapply(numeric(d), function(x) r_test)
    r_mat <- as.matrix(expand.grid(r_list))
    r_length <- dim(r_mat)[1]
    x_test[[idx]] <- apply(
      r_mat,
      1,
      funcs[[idx]]
    ) %>%
      matrix(nrow = dim(r_mat)[1], byrow = TRUE)

    r[[idx]] <- cbind(time_points[idx], r_mat)
  }

  r_full <- reduce(r, rbind)
  r_df <- data.frame(r_full)
  rnames <- paste0("r", 1:(dim(r_full)[2] - 1))
  names(r_df) <- c("time", rnames)
  x_test <- reduce(x_test, rbind)
  x_test_df <- data.frame(x_test)
  xnames <- paste0("x", 1:(dim(x_test)[2]))
  names(x_test_df) <- xnames

  full_df <- bind_cols(r_df, x_test_df)

  return(list(pme_results, unlist(msd), full_df))
}
