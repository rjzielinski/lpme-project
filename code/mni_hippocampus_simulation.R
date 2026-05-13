library(divest)
library(doFuture)
library(dplyr)
library(foreach)
library(fslr)
library(here)
library(oro.dicom)
library(oro.nifti)
library(neurobase)
library(plotly)
library(pme)
library(pracma)
library(progressr)
library(purrr)
library(readr)
library(shapes)
library(stringr)
library(tibble)
library(tidyverse)
library(lubridate)

source(here("code/functions/calculate_lpme_reconstructions.R"))
source(here("code/functions/calculate_pme_reconstructions.R"))

mni_hipp_mask <- paste0(
  fsl_dir(),
  "/data/standard/MNI152_T1_1mm_Hipp_mask_dil8.nii.gz"
) |>
  readnii()

hipp_idx <- which(mni_hipp_mask > 0, arr.ind = TRUE)
hipp_idx <- scale(hipp_idx, scale = FALSE)

lhipp_idx <- hipp_idx[hipp_idx[, 1] <= 0, ]

surface <- rep(FALSE, nrow(lhipp_idx))

for (dim in 1:3) {
  for (dim2 in 1:3) {
    if (dim != dim2) {
      dim3 <- seq(1, 3, 1)[!(1:3 %in% c(dim, dim2))]
      # at each point in a 2-dimensional profile view of structure,
      # obtain all "depth" coordinates of non-zero intensity voxels
      unique_vals <- unique(lhipp_idx[, c(dim, dim2)])

      for (i in seq_len(nrow(unique_vals))) {
        vals <- lhipp_idx[
          lhipp_idx[, dim] == unique_vals[i, 1] &
            lhipp_idx[, dim2] == unique_vals[i, 2],
          dim3
        ]
        min_vals <- vector()
        max_vals <- vector()
        # voxel is assumed to be a surface point if there are no
        # neighboring voxels in at least one direction along dim

        for (v in vals) {
          if (!((v - 1) %in% vals)) {
            min_vals <- c(min_vals, v)
          } else if (!((v + 1) %in% vals)) {
            max_vals <- c(max_vals, v)
          }
        }

        surface_vals <- c(min_vals, max_vals)
        for (v in surface_vals) {
          surface[which(
            lhipp_idx[, dim] == unique_vals[i, 1] &
              lhipp_idx[, dim2] == unique_vals[i, 2] &
              lhipp_idx[, dim3] == v,
            arr.ind = TRUE
          )] <- TRUE
        }
      }
    } else {
      next
    }
  }
}

lhipp_mat <- lhipp_idx
lhipp_surface <- lhipp_idx[surface, ]


simulate_data <- function(
  surface_data,
  duration,
  interval,
  time_noise,
  scale_noise = 0.1
) {
  time_points <- seq(interval, duration, by = interval)
  time_noise <- runif(
    n = length(time_points),
    min = -time_noise,
    max = time_noise
  )
  time_points <- c(0, time_points + time_noise)

  # scale all values to fall between -1 and 1
  for (dim_idx in seq_len(ncol(surface_data))) {
    surface_data[, dim_idx] <- surface_data[, dim_idx] -
      mean(surface_data[, dim_idx])
    surface_data[, dim_idx] <- surface_data[, dim_idx] /
      max(abs(surface_data[, dim_idx]))
  }
  time_points <- time_points / max(time_points)

  scale_vals <- rnorm(n = length(time_points), mean = 1, sd = scale_noise)
  observed_surfaces <- list()
  true_surfaces <- list()

  for (time_idx in seq_along(time_points)) {
    surface_data_scaled <- surface_data * scale_vals[time_idx]

    observed_surfaces[[time_idx]] <- cbind(
      rep(time_points[time_idx], nrow(surface_data_scaled)),
      surface_data_scaled
    )

    true_surfaces[[time_idx]] <- cbind(
      rep(time_points[time_idx], nrow(surface_data)),
      surface_data
    )
  }

  observed_surfaces <- do.call(rbind, observed_surfaces)
  true_surfaces <- do.call(rbind, true_surfaces)

  list(
    observed = observed_surfaces,
    true = true_surfaces
  )
}


fit_models <- function(sim_data) {
  observed_surfaces <- sim_data$observed
  true_surfaces <- sim_data$true

  time_points <- sort(unique(true_surfaces[, 1]))

  sim_lpme <- list()
  sim_pme <- list()
  lpme_projections <- list()
  pme_projections <- list()

  partition <- ifelse(observed_surfaces[, 4] > 0, 1, 2)
  partition_vals <- unique(partition)

  for (partition_idx in seq_along(partition_vals)) {
    observed_partition <- observed_surfaces[
      partition == partition_vals[partition_idx],
    ]

    sim_lpme[[partition_idx]] <- lpme(
      observed_partition,
      d = 2,
      min_clusters = 50,
      max_clusters = 100,
      verbose = FALSE,
      print_plots = FALSE
    )

    lpme_projections[[partition_idx]] <- calculate_lpme_reconstructions(
      sim_lpme[[partition_idx]],
      observed_partition
    )
    sim_pme[[partition_idx]] <- list()
    pme_projections[[partition_idx]] <- list()

    for (time_idx in seq_along(time_points)) {
      sim_pme[[partition_idx]][[time_idx]] <- pme(
        observed_partition[
          observed_partition[, 1] == time_points[time_idx],
          -1
        ],
        d = 2,
        lambda = exp(-20:5),
        min_clusters = 50,
        max_clusters = 100
      )

      time_pme_projections <- calculate_pme_reconstructions(
        sim_pme[[partition_idx]][[time_idx]],
        observed_partition[observed_partition[, 1] == time_points[time_idx], -1]
      )

      pme_projections[[partition_idx]][[time_idx]] <- cbind(
        rep(time_points[time_idx], nrow(time_pme_projections)),
        time_pme_projections
      )
    }

    pme_projections[[partition_idx]] <- do.call(
      rbind,
      pme_projections[[partition_idx]]
    )
  }

  lpme_projections <- do.call(rbind, lpme_projections)
  pme_projections <- do.call(rbind, pme_projections)

  observed_surfaces <- rbind(
    observed_surfaces[partition == 1, ],
    observed_surfaces[partition == 2, ]
  )

  true_surfaces <- rbind(
    true_surfaces[partition == 1, ],
    true_surfaces[partition == 2, ]
  )

  list(
    observed = observed_surfaces,
    true = true_surfaces,
    lpme_projections = lpme_projections,
    pme_projections = pme_projections,
    lpme_ests = sim_lpme,
    pme_ests = sim_pme,
    partitions = partition
  )
}


calc_procrustes <- function(estimates) {
  observed_surfaces <- estimates$observed
  true_surfaces <- estimates$true

  lpme_projections <- estimates$lpme_projections
  pme_projections <- estimates$pme_projections

  time_points <- sort(unique(true_surfaces[, 1]))

  lpme_proc_dists <- vector(mode = "numeric", length = length(time_points))
  lpme_proc_dists <- vector(mode = "numeric", length = length(time_points))
  pme_proc_dists <- vector(mode = "numeric", length = length(time_points))

  for (time_idx in seq_along(time_points)) {
    time_val <- time_points[time_idx]
    temp_true <- true_surfaces[true_surfaces[, 1] == time_val, -1]
    temp_lpme <- lpme_projections[lpme_projections[, 1] == time_val, -1]
    temp_pme <- pme_projections[pme_projections[, 1] == time_val, -1]

    lpme_proc_dists[time_idx] <- procdist(
      temp_true,
      temp_lpme,
      type = "partial"
    )

    pme_proc_dists[time_idx] <- procdist(
      temp_true,
      temp_pme,
      type = "partial"
    )
  }

  lpme_proc_dist <- sum(lpme_proc_dists)
  pme_proc_dist <- sum(pme_proc_dists)

  list(
    lpme_dist = lpme_proc_dist,
    pme_dist = pme_proc_dist
  )
}


set.seed(500)

cores <- availableCores() - 4
plan(multicore, workers = cores)

n_sims <- 1000


with_progress({
  p <- progressor(n_sims)

  lhipp_distances <- foreach(
    sim_idx = seq_len(n_sims),
    .options.future = list(seed = TRUE)
  ) %dofuture%
    {
      sim_data <- simulate_data(
        lhipp_surface,
        duration = 5,
        interval = 0.5,
        time_noise = 0.05,
        scale_noise = 0.15
      )

      lhipp_ests <- fit_models(sim_data)

      p()

      procrustes_distances <- calc_procrustes(lhipp_ests)
    }
})

saveRDS(lhipp_distances, here("output/mni_lhipp_simulation_dists.RDS"))
