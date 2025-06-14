library(alphashape3d)
library(here)
library(lubridate)
library(mirai)
library(future)
library(plotly)
library(plot3D)
library(pme)
library(pracma)
library(RColorBrewer)
library(Rfast)
library(tictoc)
library(tidyverse)

source(here("code/functions/calculate_lpme_reconstructions.R"))
source(here("code/functions/calculate_pme_reconstructions.R"))
source(here("code/functions/estimate_volume.R"))
source(here("code/functions/interior_identification.R"))
source(here("code/prinSurf_v3.R"))
source(here("code/functions/estimate_volume_interior.R"))
source(here("code/functions/preprocess_adni.R"))
source(here("code/functions/fit_adni.R"))

# options(future.globals.maxSize = 4 * 1e9)

voxel_vol <- 1.2 * 0.9375 * 0.9375

lhipp_surface <- read_csv(here("data/lhipp_surface_fsl.csv"))
rhipp_surface <- read_csv(here("data/rhipp_surface_fsl.csv"))
lthal_surface <- read_csv(here("data/lthal_surface_fsl.csv"))
rthal_surface <- read_csv(here("data/rthal_surface_fsl.csv"))

adni_info <- read_csv(here("data/adni_info_full.csv")) |>
  distinct()
adni <- read_csv(here("data/adni_info.csv")) |>
  rename(
    image_id = "Image Data ID",
    subid = "Subject",
    date = "Acq Date"
  )

lhipp_processed <- preprocess_adni(lhipp_surface, adni_info)
lhipp_surface <- lhipp_processed$surface
lhipp_centers <- lhipp_processed$centers

rhipp_processed <- preprocess_adni(rhipp_surface, adni_info)
rhipp_surface <- rhipp_processed$surface
rhipp_centers <- rhipp_processed$centers

lthal_processed <- preprocess_adni(lthal_surface, adni_info)
lthal_surface <- lthal_processed$surface
lthal_centers <- lthal_processed$centers

rthal_processed <- preprocess_adni(rthal_surface, adni_info)
rthal_surface <- rthal_processed$surface
rthal_centers <- rthal_processed$centers

write_csv(lhipp_surface, here("data/lhipp_surface_fsl_processed.csv"))
write_csv(rhipp_surface, here("data/rhipp_surface_fsl_processed.csv"))
write_csv(lthal_surface, here("data/lthal_surface_fsl_processed.csv"))
write_csv(rthal_surface, here("data/rthal_surface_fsl_processed.csv"))

patnos <- unique(lhipp_surface$subid)[1:2]

dx <- read_csv("data/DXSUM.csv")
max_dx <- dx |>
  filter(PTID %in% patnos) |>
  group_by(PTID) |>
  arrange(EXAMDATE) |>
  summarize(
    bl_dx = dplyr::first(DIAGNOSIS),
    final_dx = dplyr::last(DIAGNOSIS, na_rm = TRUE),
    max_dx = max(DIAGNOSIS, na.rm = TRUE)
  )
max_dx |>
  group_by(bl_dx) |>
  tally()


adni |>
  filter(subid %in% patnos) |>
  group_by(subid) |>
  summarize(study_group = first(Group)) |>
  ungroup() |>
  group_by(study_group) |>
  tally()

lhipp_fit <- fit_adni(
  filter(lhipp_surface, subid %in% patnos),
  filter(lhipp_centers, subid %in% patnos),
  "lhipp"
)
rhipp_fit <- fit_adni(
  filter(rhipp_surface, subid %in% patnos),
  filter(rhipp_centers, subid %in% patnos),
  "rhipp"
)
lthal_fit <- fit_adni(
  filter(lthal_surface, subid %in% patnos),
  filter(lthal_centers, subid %in% patnos),
  "lthal"
)
rthal_fit <- fit_adni(
  filter(rthal_surface, subid %in% patnos),
  filter(rthal_centers, subid %in% patnos),
  "rthal"
)

adni_files <-
  rthal_results <- list()

for (patno_idx in seq_along(patnos)) {
  patno <- patnos[patno_idx]
  patno_rthal <- rthal_surface |>
    filter(subid == patno) |>
    select(time_from_bl, x, y, z, theta, phi) |>
    mutate(
      partition = z > 0,
      partition1 = z > -0.15,
      partition2 = z < 0.15
    )

  patno_rthal_centers <- rthal_centers |>
    filter(subid == patno) |>
    group_by(date) |>
    summarize(
      max_x = max(max_x),
      max_y = max(max_y),
      max_z = max(max_z)
    )

  rthal_aug <- patno_rthal |>
    select(-partition, -partition1, -partition2) |>
    as.matrix()

  rthal_pt1 <- patno_rthal |>
    filter(partition1 == TRUE) |>
    select(time_from_bl, x, y, z) |>
    as.matrix()
  rthal_pt2 <- patno_rthal |>
    filter(partition2 == TRUE) |>
    select(time_from_bl, x, y, z) |>
    as.matrix()

  time_values <- unique(patno_rthal$time_from_bl)

  rthal_results[[patno_idx]] <- future(
    {
      rthal_pme_aug_list <- list()
      rthal_pme_aug_time_list <- list()
      rthal_pme_aug_reconstruction_list <- list()

      rthal_lpme_part <- replicate(2, NULL, simplify = FALSE)
      rthal_pme_part <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      rthal_pc_part <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      rthal_lpme_part_reconstructions <- replicate(2, NULL, simplify = FALSE)
      rthal_pme_part_reconstructions <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      rthal_pc_part_reconstructions <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      rthal_lpme_part_times <- replicate(2, NULL, simplify = FALSE)
      rthal_pme_part_times <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      rthal_pc_part_times <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      tic()
      rthal_lpme_aug <- lpme(
        rthal_aug,
        d = 2,
        verbose = TRUE,
        print_plots = FALSE
      )
      rthal_lpme_aug_time <- toc()
      rthal_lpme_aug_time <- rthal_lpme_aug_time$toc - rthal_lpme_aug_time$tic

      rthal_lpme_aug_reconstructions <- calculate_lpme_reconstructions(
        rthal_lpme_aug,
        rthal_aug
      )

      tic()
      rthal_lpme_part[[1]] <- lpme(
        rthal_pt1,
        d = 2,
        verbose = FALSE,
        print_plots = FALSE
      )
      rthal_lpme_time_pt1 <- toc()
      rthal_lpme_part_times[[1]] <- rthal_lpme_time_pt1$toc -
        rthal_lpme_time_pt1$tic
      rthal_lpme_part_reconstructions[[1]] <- calculate_lpme_reconstructions(
        rthal_lpme_part[[1]],
        rthal_pt1[rthal_pt1[, 4] > 0, ]
      )
      tic()
      rthal_lpme_part[[2]] <- lpme(
        rthal_pt2,
        d = 2,
        verbose = FALSE,
        print_plots = FALSE
      )
      rthal_lpme_time_pt2 <- toc()
      rthal_lpme_part_times[[2]] <- rthal_lpme_time_pt2$toc -
        rthal_lpme_time_pt2$tic
      rthal_lpme_part_reconstructions[[2]] <- calculate_lpme_reconstructions(
        rthal_lpme_part[[2]],
        rthal_pt2[rthal_pt2[, 4] <= 0, ]
      )

      for (time_idx in seq_along(time_values)) {
        temp_rthal <- rthal_aug[rthal_aug[, 1] == time_values[time_idx], -1]
        temp_rthal_pt1 <- rthal_pt1[rthal_pt1[, 1] == time_values[time_idx], -1]
        temp_rthal_pt2 <- rthal_pt2[rthal_pt2[, 1] == time_values[time_idx], -1]

        tic()
        rthal_pme_aug_list[[time_idx]] <- pme(temp_rthal, d = 2)
        temp_pme_aug_time <- toc()
        rthal_pme_aug_time_list[[time_idx]] <- temp_pme_aug_time$toc -
          temp_pme_aug_time$tic
        rthal_pme_aug_reconstruction_list[[
          time_idx
        ]] <- calculate_pme_reconstructions(
          rthal_pme_aug_list[[time_idx]],
          temp_rthal
        )
        rthal_pme_aug_reconstruction_list[[time_idx]] <- cbind(
          time_values[time_idx],
          rthal_pme_aug_reconstruction_list[[time_idx]]
        )

        tic()
        rthal_pme_part[[1]][[time_idx]] <- pme(temp_rthal_pt1, d = 2)
        temp_pme_time_pt1 <- toc()
        rthal_pme_part_times[[1]][[time_idx]] <- temp_pme_time_pt1$toc -
          temp_pme_time_pt1$tic
        rthal_pme_part_reconstructions[[1]][[
          time_idx
        ]] <- calculate_pme_reconstructions(
          rthal_pme_part[[1]][[time_idx]],
          temp_rthal_pt1[temp_rthal_pt1[, 3] > 0, ]
        )
        tic()
        rthal_pme_part[[2]][[time_idx]] <- pme(temp_rthal_pt2, d = 2)
        temp_pme_time_pt2 <- toc()
        rthal_pme_part_times[[2]][[time_idx]] <- temp_pme_time_pt2$toc -
          temp_pme_time_pt2$tic
        rthal_pme_part_reconstructions[[2]][[
          time_idx
        ]] <- calculate_pme_reconstructions(
          rthal_pme_part[[2]][[time_idx]],
          temp_rthal_pt2[temp_rthal_pt2[, 3] <= 0, ]
        )

        rthal_pme_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          rthal_pme_part_reconstructions[[1]][[time_idx]]
        )
        rthal_pme_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          rthal_pme_part_reconstructions[[2]][[time_idx]]
        )

        tic()
        principal_surface_part1 <- prinSurf(temp_rthal_pt1)
        temp_pc_time_pt1 <- toc()
        rthal_pc_part_times[[1]][[time_idx]] <- temp_pc_time_pt1$toc -
          temp_pc_time_pt1$tic
        surface_mse_part1 <- map(
          seq_along(principal_surface_part1),
          ~ principal_surface_part1[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part1 <- which.min(surface_mse_part1)

        rthal_pc_part[[1]][[time_idx]] <- principal_surface_part1[[
          opt_surface_part1 + 2
        ]]
        rthal_pc_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part1[[opt_surface_part1 + 2]]$PS
        )
        rthal_pc_part_reconstructions[[1]][[
          time_idx
        ]] <- rthal_pc_part_reconstructions[[1]][[time_idx]][
          temp_rthal_pt1[, 3] > 0,
        ]

        tic()
        principal_surface_part2 <- prinSurf(temp_rthal_pt2)
        temp_pc_time_pt2 <- toc()
        rthal_pc_part_times[[2]][[time_idx]] <- temp_pc_time_pt2$toc -
          temp_pc_time_pt2$tic
        surface_mse_part2 <- map(
          seq_along(principal_surface_part2),
          ~ principal_surface_part2[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part2 <- which.min(surface_mse_part2)

        rthal_pc_part[[2]][[time_idx]] <- principal_surface_part2[[
          opt_surface_part2 + 2
        ]]
        rthal_pc_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part2[[opt_surface_part2 + 2]]$PS
        )
        rthal_pc_part_reconstructions[[2]][[
          time_idx
        ]] <- rthal_pc_part_reconstructions[[2]][[time_idx]][
          temp_rthal_pt2[, 3] <= 0,
        ]
      }

      rthal_lpme_part_reconstructions <- reduce(
        rthal_lpme_part_reconstructions,
        rbind
      )

      rthal_lpme_part_times <- reduce(
        rthal_lpme_part_times,
        sum
      )

      rthal_pme_aug_reconstructions <- reduce(
        rthal_pme_aug_reconstruction_list,
        rbind
      )

      rthal_pme_aug_time_list <- reduce(
        rthal_pme_aug_time_list,
        sum
      )

      rthal_pme_part_reconstructions <- map(
        seq_along(rthal_pme_part_reconstructions),
        ~ reduce(rthal_pme_part_reconstructions[[.x]], rbind)
      )
      rthal_pme_part_reconstructions <- reduce(
        rthal_pme_part_reconstructions,
        rbind
      )

      rthal_pme_part_times <- map(
        seq_along(rthal_pme_part_times),
        ~ reduce(rthal_pme_part_times[[.x]], sum)
      )
      rthal_pme_part_times <- reduce(
        rthal_pme_part_times,
        sum
      )

      rthal_pc_part_reconstructions <- map(
        seq_along(rthal_pc_part_reconstructions),
        ~ reduce(rthal_pc_part_reconstructions[[.x]], rbind)
      )
      rthal_pc_part_reconstructions <- reduce(
        rthal_pc_part_reconstructions,
        rbind
      )

      rthal_pc_part_times <- map(
        seq_along(rthal_pc_part_times),
        ~ reduce(rthal_pc_part_times[[.x]], sum)
      )
      rthal_pc_part_times <- reduce(
        rthal_pc_part_times,
        sum
      )

      rthal_lpme_part_volumes <- estimate_volume_interior_lpme(
        rthal_lpme_part,
        list(rthal_pt1, rthal_pt2),
        time_values,
        n_points = 10000,
        data_max = patno_rthal_centers,
        limit_scaler = 0.05,
        partition_index = 3
      )
      rthal_lpme_part_volumes <- rthal_lpme_part_volumes$volumes

      rthal_pme_part_volumes <- estimate_volume_interior_pme(
        rthal_pme_part,
        list(rthal_pt1, rthal_pt2),
        time_values,
        n_points = 10000,
        data_max = patno_rthal_centers,
        limit_scaler = 0.05,
        partition_index = 3
      )
      rthal_pme_part_volumes <- rthal_pme_part_volumes$volumes

      rthal_lpme_aug_volumes <- vector()
      rthal_pme_aug_volumes <- vector()
      rthal_pc_part_volumes <- vector()

      rthal_lpme_aug_interior <- list()
      rthal_pme_aug_interior <- list()
      rthal_pc_part_interior <- list()

      rthal_lpme_part_volumes_mesh <- vector()
      rthal_pme_part_volumes_mesh <- vector()

      for (time_idx in seq_along(time_values)) {
        x_scale <- patno_rthal_centers$max_x[time_idx]
        y_scale <- patno_rthal_centers$max_y[time_idx]
        z_scale <- patno_rthal_centers$max_z[time_idx]
        temp_lpme_reconstructions <- rthal_lpme_aug_reconstructions[
          rthal_lpme_aug_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pme_reconstructions <- rthal_pme_aug_reconstructions[
          rthal_pme_aug_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pc_reconstructions <- rthal_pc_part_reconstructions[
          rthal_pc_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]

        temp_lpme_part_reconstructions <- rthal_lpme_part_reconstructions[
          rthal_lpme_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pme_part_reconstructions <- rthal_pme_part_reconstructions[
          rthal_pme_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]

        temp_lpme_ashape <- ashape3d(
          temp_lpme_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        rthal_lpme_aug_volumes[time_idx] <- volume_ashape3d(
          temp_lpme_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pme_ashape <- ashape3d(
          temp_pme_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        rthal_pme_aug_volumes[time_idx] <- volume_ashape3d(
          temp_pme_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pc_ashape <- ashape3d(
          temp_pc_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        rthal_pc_part_volumes[time_idx] <- volume_ashape3d(
          temp_pc_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_lpme_part_ashape <- ashape3d(
          temp_lpme_part_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        rthal_lpme_part_volumes_mesh[time_idx] <- volume_ashape3d(
          temp_lpme_part_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pme_part_ashape <- ashape3d(
          temp_pme_part_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        rthal_pme_part_volumes_mesh[time_idx] <- volume_ashape3d(
          temp_pme_part_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)
      }

      rthal_lpme_aug_out <- list(
        lpme = rthal_lpme_aug,
        reconstructions = rthal_lpme_aug_reconstructions,
        volumes = rthal_lpme_aug_volumes,
        fit_time = rthal_lpme_aug_time
      )
      rthal_pme_aug_out <- list(
        pme = rthal_pme_aug_list,
        reconstructions = rthal_pme_aug_reconstructions,
        volumes = rthal_pme_aug_volumes,
        fit_time = rthal_pme_aug_time_list
      )

      rthal_lpme_part_out <- list(
        lpme = rthal_lpme_part,
        reconstructions = rthal_lpme_part_reconstructions,
        volumes = rthal_lpme_part_volumes,
        volumes_mesh = rthal_lpme_part_volumes_mesh,
        fit_time = rthal_lpme_part_times
      )

      rthal_pme_part_out <- list(
        pme = rthal_pme_part,
        reconstructions = rthal_pme_part_reconstructions,
        volumes = rthal_pme_part_volumes,
        volumes_mesh = rthal_pme_part_volumes_mesh,
        fit_time = rthal_pme_part_times
      )

      rthal_pc_part_out <- list(
        pc = rthal_pc_part,
        reconstructions = rthal_pc_part_reconstructions,
        volumes = rthal_pc_part_volumes,
        fit_time = rthal_pc_part_times
      )

      rthal <- list(
        lpme_aug = rthal_lpme_aug_out,
        pme_aug = rthal_pme_aug_out,
        lpme_part = rthal_lpme_part_out,
        pme_part = rthal_pme_part_out,
        pc_part = rthal_pc_part_out
      )

      if (!dir.exists(here(paste0("output/adni/", patno)))) {
        dir.create(here(paste0("output/adni/", patno)), recursive = TRUE)
      }

      saveRDS(
        rthal,
        here(paste0("output/adni/", patno, "/rthal_results.RDS"))
      )

      rthal_lpme_aug_out <- NULL
      rthal_pme_aug_out <- NULL
      rthal_lpme_part_out <- NULL
      rthal_pme_part_out <- NULL
      rthal_pc_part_out <- NULL

      gc()

      rthal_out <- tibble(
        patno = patno,
        rthal_time_values = time_values,
        rthal_lpme_aug_vols = rthal_lpme_aug_volumes,
        rthal_pme_aug_vols = rthal_pme_aug_volumes,
        rthal_lpme_part_vols = rthal_lpme_part_volumes,
        rthal_pme_part_vols = rthal_pme_part_volumes,
        rthal_pc_part_vols = rthal_pc_part_volumes,
        rthal_lpme_part_vols_mesh = rthal_lpme_part_volumes_mesh,
        rthal_pme_part_vols_mesh = rthal_pme_part_volumes_mesh,
        rthal_lpme_aug_time = rthal_lpme_aug_time,
        rthal_pme_aug_time = rthal_pme_aug_time_list,
        rthal_lpme_part_time = rthal_lpme_part_times,
        rthal_pme_part_time = rthal_pme_part_times,
        rthal_pc_part_time = rthal_pc_part_times
      )
      rthal_out
    },
    seed = TRUE
  )
}

rthal_results_out <- purrr::map(rthal_results, value)

# rthal_df_list <- list()

rthal_df <- bind_rows(rthal_results_out)

adni_volumes <- bind_cols(
  lhipp_df,
  rhipp_df,
  lthal_df,
  rthal_df
)

if (!dir.exists(here("output/adni/", patno))) {
  dir.create(here("output/adni/"), recursive = TRUE)
}

write_csv(adni_volumes, here("output/adni/adni_volumes.csv"))
