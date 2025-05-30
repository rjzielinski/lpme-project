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

options(future.globals.maxSize = 4 * 1e9)

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

lhipp_surface <- lhipp_surface |>
  left_join(adni_info, by = join_by(scan_id == image_id))
rhipp_surface <- rhipp_surface |>
  left_join(adni_info, by = join_by(scan_id == image_id))
lthal_surface <- lthal_surface |>
  left_join(adni_info, by = join_by(scan_id == image_id))
rthal_surface <- rthal_surface |>
  left_join(adni_info, by = join_by(scan_id == image_id))

voxel_vol <- 1.2 * 0.9375 * 0.9375


lhipp_surface <- lhipp_surface |>
  mutate(
    date = decimal_date(date)
  )

rhipp_surface <- rhipp_surface |>
  mutate(
    date = decimal_date(date)
  )

lthal_surface <- lthal_surface |>
  mutate(
    date = decimal_date(date)
  )

rthal_surface <- rthal_surface |>
  mutate(
    date = decimal_date(date)
  )


lhipp_centers <- lhipp_surface |>
  group_by(subid, date, scan_id) |>
  summarise(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x - mean_x)),
    max_y = max(abs(y - mean_y)),
    max_z = max(abs(z - mean_z))
  ) |>
  ungroup()

rhipp_centers <- rhipp_surface |>
  group_by(subid, date, scan_id) |>
  summarise(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x - mean_x)),
    max_y = max(abs(y - mean_y)),
    max_z = max(abs(z - mean_z))
  ) |>
  ungroup()

lthal_centers <- lthal_surface |>
  group_by(subid, date, scan_id) |>
  summarise(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x - mean_x)),
    max_y = max(abs(y - mean_y)),
    max_z = max(abs(z - mean_z))
  ) |>
  ungroup()

rthal_centers <- rthal_surface |>
  group_by(subid, date, scan_id) |>
  summarise(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x - mean_x)),
    max_y = max(abs(y - mean_y)),
    max_z = max(abs(z - mean_z))
  ) |>
  ungroup()


lhipp_bl <- lhipp_surface |>
  group_by(subid) |>
  arrange(date) |>
  summarize(
    date_bl = first(date),
    date_max = max(date)
  ) |>
  mutate(duration = date_max - date_bl)

rhipp_bl <- rhipp_surface |>
  group_by(subid) |>
  arrange(date) |>
  summarize(
    date_bl = first(date),
    date_max = max(date)
  ) |>
  mutate(duration = date_max - date_bl)

lthal_bl <- lthal_surface |>
  group_by(subid) |>
  arrange(date) |>
  summarize(
    date_bl = first(date),
    date_max = max(date)
  ) |>
  mutate(duration = date_max - date_bl)

rthal_bl <- rthal_surface |>
  group_by(subid) |>
  arrange(date) |>
  summarize(
    date_bl = first(date),
    date_max = max(date)
  ) |>
  mutate(duration = date_max - date_bl)


lhipp_surface <- lhipp_surface |>
  left_join(lhipp_centers, by = c("subid", "date", "scan_id")) |>
  left_join(lhipp_bl, by = "subid") |>
  mutate(
    x = (x - mean_x) / max_x,
    y = (y - mean_y) / max_y,
    z = (z - mean_z) / max_z,
    time_from_bl = (date - date_bl) / duration
  ) |>
  filter(duration > 2)

rhipp_surface <- rhipp_surface |>
  left_join(rhipp_centers, by = c("subid", "date", "scan_id")) |>
  left_join(rhipp_bl, by = "subid") |>
  mutate(
    x = (x - mean_x) / max_x,
    y = (y - mean_y) / max_y,
    z = (z - mean_z) / max_z,
    time_from_bl = (date - date_bl) / duration
  ) |>
  filter(duration > 2)

lthal_surface <- lthal_surface |>
  left_join(lthal_centers, by = c("subid", "date", "scan_id")) |>
  left_join(lthal_bl, by = "subid") |>
  mutate(
    x = (x - mean_x) / max_x,
    y = (y - mean_y) / max_y,
    z = (z - mean_z) / max_z,
    time_from_bl = (date - date_bl) / duration
  ) |>
  filter(duration > 2)

rthal_surface <- rthal_surface |>
  left_join(rthal_centers, by = c("subid", "date", "scan_id")) |>
  left_join(rthal_bl, by = "subid") |>
  mutate(
    x = (x - mean_x) / max_x,
    y = (y - mean_y) / max_y,
    z = (z - mean_z) / max_z,
    time_from_bl = (date - date_bl) / duration
  ) |>
  filter(duration > 2)


lhipp_surface_spherical <- lhipp_surface |>
  select(x, y, z) |>
  as.matrix() |>
  cart2sph() |>
  as_tibble()

rhipp_surface_spherical <- rhipp_surface |>
  select(x, y, z) |>
  as.matrix() |>
  cart2sph() |>
  as_tibble()

lthal_surface_spherical <- lthal_surface |>
  select(x, y, z) |>
  as.matrix() |>
  cart2sph() |>
  as_tibble()

rthal_surface_spherical <- rthal_surface |>
  select(x, y, z) |>
  as.matrix() |>
  cart2sph() |>
  as_tibble()


lhipp_surface <- bind_cols(lhipp_surface, lhipp_surface_spherical)
rhipp_surface <- bind_cols(rhipp_surface, rhipp_surface_spherical)
lthal_surface <- bind_cols(lthal_surface, lthal_surface_spherical)
rthal_surface <- bind_cols(rthal_surface, rthal_surface_spherical)


write_csv(lhipp_surface, here("data/lhipp_surface_fsl_processed.csv"))
write_csv(rhipp_surface, here("data/rhipp_surface_fsl_processed.csv"))
write_csv(lthal_surface, here("data/lthal_surface_fsl_processed.csv"))
write_csv(rthal_surface, here("data/rthal_surface_fsl_processed.csv"))

patnos <- unique(lhipp_surface$subid)
cores <- parallel::detectCores()
# plan(multisession, workers = cores)
daemons(cores)

lhipp_results <- list()

for (patno_idx in seq_along(patnos)) {
  lhipp_results[[patno_idx]] <- mirai(
    {
      patno <- patnos[patno_idx]
      patno_lhipp <- lhipp_surface |>
        filter(subid == patno) |>
        select(time_from_bl, x, y, z, theta, phi) |>
        mutate(
          partition = z > 0,
          partition1 = z > -0.15,
          partition2 = z < 0.15
        )

      patno_lhipp_centers <- lhipp_centers |>
        filter(subid == patno) |>
        group_by(date) |>
        summarize(
          max_x = max(max_x),
          max_y = max(max_y),
          max_z = max(max_z)
        )

      lhipp_aug <- patno_lhipp |>
        select(-partition, -partition1, -partition2) |>
        as.matrix()

      lhipp_pt1 <- patno_lhipp |>
        filter(partition1 == TRUE) |>
        select(time_from_bl, x, y, z) |>
        as.matrix()
      lhipp_pt2 <- patno_lhipp |>
        filter(partition2 == TRUE) |>
        select(time_from_bl, x, y, z) |>
        as.matrix()

      time_values <- unique(patno_lhipp$time_from_bl)

      lhipp_pme_aug_list <- list()
      lhipp_pme_aug_time_list <- list()
      lhipp_pme_aug_reconstruction_list <- list()

      lhipp_lpme_part <- replicate(2, NULL, simplify = FALSE)
      lhipp_pme_part <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      lhipp_pc_part <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      lhipp_lpme_part_reconstructions <- replicate(2, NULL, simplify = FALSE)
      lhipp_pme_part_reconstructions <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      lhipp_pc_part_reconstructions <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      lhipp_lpme_part_times <- replicate(2, NULL, simplify = FALSE)
      lhipp_pme_part_times <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      lhipp_pc_part_times <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      tic()
      lhipp_lpme_aug <- lpme(
        lhipp_aug,
        d = 2,
        verbose = TRUE,
        print_plots = FALSE
      )
      lhipp_lpme_aug_time <- toc()
      lhipp_lpme_aug_time <- lhipp_lpme_aug_time$toc - lhipp_lpme_aug_time$tic

      lhipp_lpme_aug_reconstructions <- calculate_lpme_reconstructions(
        lhipp_lpme_aug,
        lhipp_aug
      )

      tic()
      lhipp_lpme_part[[1]] <- lpme(
        lhipp_pt1,
        d = 2,
        verbose = FALSE,
        print_plots = FALSE
      )
      lhipp_lpme_time_pt1 <- toc()
      lhipp_lpme_part_times[[1]] <- lhipp_lpme_time_pt1$toc -
        lhipp_lpme_time_pt1$tic
      lhipp_lpme_part_reconstructions[[1]] <- calculate_lpme_reconstructions(
        lhipp_lpme_part[[1]],
        lhipp_pt1[lhipp_pt1[, 4] > 0, ]
      )
      tic()
      lhipp_lpme_part[[2]] <- lpme(
        lhipp_pt2,
        d = 2,
        verbose = FALSE,
        print_plots = FALSE
      )
      lhipp_lpme_time_pt2 <- toc()
      lhipp_lpme_part_times[[2]] <- lhipp_lpme_time_pt2$toc -
        lhipp_lpme_time_pt2$tic
      lhipp_lpme_part_reconstructions[[2]] <- calculate_lpme_reconstructions(
        lhipp_lpme_part[[2]],
        lhipp_pt2[lhipp_pt2[, 4] <= 0, ]
      )

      for (time_idx in seq_along(time_values)) {
        temp_lhipp <- lhipp_aug[lhipp_aug[, 1] == time_values[time_idx], -1]
        temp_lhipp_pt1 <- lhipp_pt1[lhipp_pt1[, 1] == time_values[time_idx], -1]
        temp_lhipp_pt2 <- lhipp_pt2[lhipp_pt2[, 1] == time_values[time_idx], -1]

        tic()
        lhipp_pme_aug_list[[time_idx]] <- pme(temp_lhipp, d = 2)
        temp_pme_aug_time <- toc()
        lhipp_pme_aug_time_list[[time_idx]] <- temp_pme_aug_time$toc -
          temp_pme_aug_time$tic
        lhipp_pme_aug_reconstruction_list[[
          time_idx
        ]] <- calculate_pme_reconstructions(
          lhipp_pme_aug_list[[time_idx]],
          temp_lhipp
        )
        lhipp_pme_aug_reconstruction_list[[time_idx]] <- cbind(
          time_values[time_idx],
          lhipp_pme_aug_reconstruction_list[[time_idx]]
        )

        tic()
        lhipp_pme_part[[1]][[time_idx]] <- pme(temp_lhipp_pt1, d = 2)
        temp_pme_time_pt1 <- toc()
        lhipp_pme_part_times[[1]][[time_idx]] <- temp_pme_time_pt1$toc -
          temp_pme_time_pt1$tic
        lhipp_pme_part_reconstructions[[1]][[
          time_idx
        ]] <- calculate_pme_reconstructions(
          lhipp_pme_part[[1]][[time_idx]],
          temp_lhipp_pt1[temp_lhipp_pt1[, 3] > 0, ]
        )
        tic()
        lhipp_pme_part[[2]][[time_idx]] <- pme(temp_lhipp_pt2, d = 2)
        temp_pme_time_pt2 <- toc()
        lhipp_pme_part_times[[2]][[time_idx]] <- temp_pme_time_pt2$toc -
          temp_pme_time_pt2$tic
        lhipp_pme_part_reconstructions[[2]][[
          time_idx
        ]] <- calculate_pme_reconstructions(
          lhipp_pme_part[[2]][[time_idx]],
          temp_lhipp_pt2[temp_lhipp_pt2[, 3] <= 0, ]
        )

        lhipp_pme_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          lhipp_pme_part_reconstructions[[1]][[time_idx]]
        )
        lhipp_pme_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          lhipp_pme_part_reconstructions[[2]][[time_idx]]
        )

        tic()
        principal_surface_part1 <- prinSurf(temp_lhipp_pt1)
        temp_pc_time_pt1 <- toc()
        lhipp_pc_part_times[[1]][[time_idx]] <- temp_pc_time_pt1$toc -
          temp_pc_time_pt1$tic
        surface_mse_part1 <- map(
          seq_along(principal_surface_part1),
          ~ principal_surface_part1[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part1 <- which.min(surface_mse_part1)

        lhipp_pc_part[[1]][[time_idx]] <- principal_surface_part1[[
          opt_surface_part1 + 2
        ]]
        lhipp_pc_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part1[[opt_surface_part1 + 2]]$PS
        )
        lhipp_pc_part_reconstructions[[1]][[
          time_idx
        ]] <- lhipp_pc_part_reconstructions[[1]][[time_idx]][
          temp_lhipp_pt1[, 3] > 0,
        ]

        tic()
        principal_surface_part2 <- prinSurf(temp_lhipp_pt2)
        temp_pc_time_pt2 <- toc()
        lhipp_pc_part_times[[2]][[time_idx]] <- temp_pc_time_pt2$toc -
          temp_pc_time_pt2$tic
        surface_mse_part2 <- map(
          seq_along(principal_surface_part2),
          ~ principal_surface_part2[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part2 <- which.min(surface_mse_part2)

        lhipp_pc_part[[2]][[time_idx]] <- principal_surface_part2[[
          opt_surface_part2 + 2
        ]]
        lhipp_pc_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part2[[opt_surface_part2 + 2]]$PS
        )
        lhipp_pc_part_reconstructions[[2]][[
          time_idx
        ]] <- lhipp_pc_part_reconstructions[[2]][[time_idx]][
          temp_lhipp_pt2[, 3] <= 0,
        ]
      }

      lhipp_lpme_part_reconstructions <- reduce(
        lhipp_lpme_part_reconstructions,
        rbind
      )

      lhipp_lpme_part_times <- reduce(
        lhipp_lpme_part_times,
        sum
      )

      lhipp_pme_aug_reconstructions <- reduce(
        lhipp_pme_aug_reconstruction_list,
        rbind
      )

      lhipp_pme_aug_time_list <- reduce(
        lhipp_pme_aug_time_list,
        sum
      )

      lhipp_pme_part_reconstructions <- map(
        seq_along(lhipp_pme_part_reconstructions),
        ~ reduce(lhipp_pme_part_reconstructions[[.x]], rbind)
      )
      lhipp_pme_part_reconstructions <- reduce(
        lhipp_pme_part_reconstructions,
        rbind
      )

      lhipp_pme_part_times <- map(
        seq_along(lhipp_pme_part_times),
        ~ reduce(lhipp_pme_part_times[[.x]], sum)
      )
      lhipp_pme_part_times <- reduce(
        lhipp_pme_part_times,
        sum
      )

      lhipp_pc_part_reconstructions <- map(
        seq_along(lhipp_pc_part_reconstructions),
        ~ reduce(lhipp_pc_part_reconstructions[[.x]], rbind)
      )
      lhipp_pc_part_reconstructions <- reduce(
        lhipp_pc_part_reconstructions,
        rbind
      )

      lhipp_pc_part_times <- map(
        seq_along(lhipp_pc_part_times),
        ~ reduce(lhipp_pc_part_times[[.x]], sum)
      )
      lhipp_pc_part_times <- reduce(
        lhipp_pc_part_times,
        sum
      )

      lhipp_lpme_part_volumes <- estimate_volume_interior_lpme(
        lhipp_lpme_part,
        list(lhipp_pt1, lhipp_pt2),
        time_values,
        n_points = 10000,
        data_max = patno_lhipp_centers,
        limit_scaler = 0.05
      )
      lhipp_lpme_part_volumes <- lhipp_lpme_part_volumes$volumes

      lhipp_pme_part_volumes <- estimate_volume_interior_pme(
        lhipp_pme_part,
        list(lhipp_pt1, lhipp_pt2),
        time_values,
        n_points = 10000,
        data_max = patno_lhipp_centers,
        limit_scaler = 0.05
      )
      lhipp_pme_part_volumes <- lhipp_pme_part_volumes$volumes

      lhipp_lpme_aug_volumes <- vector()
      lhipp_pme_aug_volumes <- vector()
      lhipp_pc_part_volumes <- vector()

      lhipp_lpme_aug_interior <- list()
      lhipp_pme_aug_interior <- list()
      lhipp_pc_part_interior <- list()

      lhipp_lpme_part_volumes_mesh <- vector()
      lhipp_pme_part_volumes_mesh <- vector()

      for (time_idx in seq_along(time_values)) {
        x_scale <- patno_lhipp_centers$max_x[time_idx]
        y_scale <- patno_lhipp_centers$max_y[time_idx]
        z_scale <- patno_lhipp_centers$max_z[time_idx]
        temp_lpme_reconstructions <- lhipp_lpme_aug_reconstructions[
          lhipp_lpme_aug_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pme_reconstructions <- lhipp_pme_aug_reconstructions[
          lhipp_pme_aug_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pc_reconstructions <- lhipp_pc_part_reconstructions[
          lhipp_pc_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]

        temp_lpme_part_reconstructions <- lhipp_lpme_part_reconstructions[
          lhipp_lpme_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pme_part_reconstructions <- lhipp_pme_part_reconstructions[
          lhipp_pme_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]

        temp_lpme_ashape <- ashape3d(
          temp_lpme_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        lhipp_lpme_aug_volumes[time_idx] <- volume_ashape3d(
          temp_lpme_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pme_ashape <- ashape3d(
          temp_pme_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        lhipp_pme_aug_volumes[time_idx] <- volume_ashape3d(
          temp_pme_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pc_ashape <- ashape3d(
          temp_pc_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        lhipp_pc_part_volumes[time_idx] <- volume_ashape3d(
          temp_pc_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_lpme_part_ashape <- ashape3d(
          temp_lpme_part_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        lhipp_lpme_part_volumes_mesh[time_idx] <- volume_ashape3d(
          temp_lpme_part_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pme_part_ashape <- ashape3d(
          temp_pme_part_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        lhipp_pme_part_volumes_mesh[time_idx] <- volume_ashape3d(
          temp_pme_part_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)
      }

      lhipp_lpme_aug_out <- list(
        lpme = lhipp_lpme_aug,
        reconstructions = lhipp_lpme_aug_reconstructions,
        volumes = lhipp_lpme_aug_volumes,
        fit_time = lhipp_lpme_aug_time
      )
      lhipp_pme_aug_out <- list(
        pme = lhipp_pme_aug_list,
        reconstructions = lhipp_pme_aug_reconstructions,
        volumes = lhipp_pme_aug_volumes,
        fit_time = lhipp_pme_aug_time_list
      )

      lhipp_lpme_part_out <- list(
        lpme = lhipp_lpme_part,
        reconstructions = lhipp_lpme_part_reconstructions,
        volumes = lhipp_lpme_part_volumes,
        volumes_mesh = lhipp_lpme_part_volumes_mesh,
        fit_time = lhipp_lpme_part_times
      )

      lhipp_pme_part_out <- list(
        pme = lhipp_pme_part,
        reconstructions = lhipp_pme_part_reconstructions,
        volumes = lhipp_pme_part_volumes,
        volumes_mesh = lhipp_pme_part_volumes_mesh,
        fit_time = lhipp_pme_part_times
      )

      lhipp_pc_part_out <- list(
        pc = lhipp_pc_part,
        reconstructions = lhipp_pc_part_reconstructions,
        volumes = lhipp_pc_part_volumes,
        fit_time = lhipp_pc_part_times
      )

      models <- list(
        lpme_aug = lhipp_lpme_aug_out,
        pme_aug = lhipp_pme_aug_out,
        lpme_part = lhipp_lpme_part_out,
        pme_part = lhipp_pme_part_out,
        pc_part = lhipp_pc_part_out
      )

      if (!dir.exists(here(paste0("output/adni/", patno)))) {
        dir.create(here(paste0("output/adni/", patno)), recursive = TRUE)
      }

      saveRDS(
        models,
        here(paste0("output/adni/", patno, "/lhipp_results.RDS"))
      )

      models
    },
    seed = TRUE,
    patnos = patnos,
    patno_idx = patno_idx,
    lhipp_surface = lhipp_surface,
    lhipp_centers = lhipp_centers
  )
}

# lhipp_results_out <- purrr::map(lhipp_results, value)
lhipp_results_out <- purrr::map(lhipp_results, ~ .x[])


lhipp_df_list <- list()

for (patno_idx in seq_along(patnos)) {
  patno <- patnos[patno_idx]
  lhipp <- lhipp_results_out[[patno_idx]]

  lhipp_lpme_aug <- lhipp$lpme_aug
  lhipp_pme_aug <- lhipp$pme_aug
  lhipp_lpme_part <- lhipp$lpme_part
  lhipp_pme_part <- lhipp$pme_part
  lhipp_pc_part <- lhipp$pc_part

  lhipp_time_values <- unique(lhipp_lpme_aug$reconstructions[, 1])
  lhipp_lpme_aug_vols <- lhipp_lpme_aug$volumes$volumes
  lhipp_lpme_aug_time <- lhipp_lpme_aug$fit_time

  lhipp_pme_aug_vols <- lhipp_pme_aug$volumes$volumes
  lhipp_pme_aug_time <- lhipp_pme_aug$fit_time

  lhipp_lpme_part_vols <- lhipp_lpme_part$volumes$volumes
  lhipp_lpme_part_vols_mesh <- lhipp_lpme_part$volumes_mesh
  lhipp_lpme_part_time <- lhipp_lpme_part$fit_time

  lhipp_pme_part_vols <- lhipp_pme_part$volumes$volumes
  lhipp_pme_part_vols_mesh <- lhipp_pme_part$volumes_mesh
  lhipp_pme_part_time <- lhipp_pme_part$fit_time

  lhipp_pc_part_vols <- lhipp_pc_part$volumes$volumes
  lhipp_pc_part_time <- lhipp_pc_part$fit_time

  lhipp_df_list[[patno_idx]] <- tibble(
    patno = patno,
    lhipp_time_values = lhipp_time_values,
    lhipp_lpme_aug_vols = lhipp_lpme_aug_vols,
    lhipp_pme_aug_vols = lhipp_pme_aug_vols,
    lhipp_lpme_part_vols = lhipp_lpme_part_vols,
    lhipp_pme_part_vols = lhipp_pme_part_vols,
    lhipp_pc_part_vols = lhipp_pc_part_vols,
    lhipp_lpme_part_vols_mesh = lhipp_lpme_part_vols_mesh,
    lhipp_pme_part_vols_mesh = lhipp_pme_part_vols_mesh,
    lhipp_lpme_aug_time = lhipp_lpme_aug_time,
    lhipp_pme_aug_time = lhipp_pme_aug_time,
    lhipp_lpme_part_time = lhipp_lpme_part_time,
    lhipp_pme_part_time = lhipp_pme_part_time,
    lhipp_pc_part_time = lhipp_pc_part_time
  )
}

lhipp_df <- bind_rows(lhipp_df_list)

rhipp_results <- list()

for (patno_idx in seq_along(patnos)) {
  rhipp_results[[patno_idx]] <- mirai(
    {
      patno <- patnos[patno_idx]
      patno_rhipp <- rhipp_surface |>
        filter(subid == patno) |>
        select(time_from_bl, x, y, z, theta, phi) |>
        mutate(
          partition = z > 0,
          partition1 = z > -0.15,
          partition2 = z < 0.15
        )

      patno_rhipp_centers <- rhipp_centers |>
        filter(subid == patno) |>
        group_by(date) |>
        summarize(
          max_x = max(max_x),
          max_y = max(max_y),
          max_z = max(max_z)
        )

      rhipp_aug <- patno_rhipp |>
        select(-partition, -partition1, -partition2) |>
        as.matrix()

      rhipp_pt1 <- patno_rhipp |>
        filter(partition1 == TRUE) |>
        select(time_from_bl, x, y, z) |>
        as.matrix()
      rhipp_pt2 <- patno_rhipp |>
        filter(partition2 == TRUE) |>
        select(time_from_bl, x, y, z) |>
        as.matrix()

      time_values <- unique(patno_rhipp$time_from_bl)

      rhipp_pme_aug_list <- list()
      rhipp_pme_aug_time_list <- list()
      rhipp_pme_aug_reconstruction_list <- list()

      rhipp_lpme_part <- replicate(2, NULL, simplify = FALSE)
      rhipp_pme_part <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      rhipp_pc_part <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      rhipp_lpme_part_reconstructions <- replicate(2, NULL, simplify = FALSE)
      rhipp_pme_part_reconstructions <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      rhipp_pc_part_reconstructions <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      rhipp_lpme_part_times <- replicate(2, NULL, simplify = FALSE)
      rhipp_pme_part_times <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      rhipp_pc_part_times <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      tic()
      rhipp_lpme_aug <- lpme(
        rhipp_aug,
        d = 2,
        verbose = TRUE,
        print_plots = FALSE
      )
      rhipp_lpme_aug_time <- toc()
      rhipp_lpme_aug_time <- rhipp_lpme_aug_time$toc - rhipp_lpme_aug_time$tic

      rhipp_lpme_aug_reconstructions <- calculate_lpme_reconstructions(
        rhipp_lpme_aug,
        rhipp_aug
      )

      tic()
      rhipp_lpme_part[[1]] <- lpme(
        rhipp_pt1,
        d = 2,
        verbose = FALSE,
        print_plots = FALSE
      )
      rhipp_lpme_time_pt1 <- toc()
      rhipp_lpme_part_times[[1]] <- rhipp_lpme_time_pt1$toc -
        rhipp_lpme_time_pt1$tic
      rhipp_lpme_part_reconstructions[[1]] <- calculate_lpme_reconstructions(
        rhipp_lpme_part[[1]],
        rhipp_pt1[rhipp_pt1[, 4] > 0, ]
      )
      tic()
      rhipp_lpme_part[[2]] <- lpme(
        rhipp_pt2,
        d = 2,
        verbose = FALSE,
        print_plots = FALSE
      )
      rhipp_lpme_time_pt2 <- toc()
      rhipp_lpme_part_times[[2]] <- rhipp_lpme_time_pt2$toc -
        rhipp_lpme_time_pt2$tic
      rhipp_lpme_part_reconstructions[[2]] <- calculate_lpme_reconstructions(
        rhipp_lpme_part[[2]],
        rhipp_pt2[rhipp_pt2[, 4] <= 0, ]
      )

      for (time_idx in seq_along(time_values)) {
        temp_rhipp <- rhipp_aug[rhipp_aug[, 1] == time_values[time_idx], -1]
        temp_rhipp_pt1 <- rhipp_pt1[rhipp_pt1[, 1] == time_values[time_idx], -1]
        temp_rhipp_pt2 <- rhipp_pt2[rhipp_pt2[, 1] == time_values[time_idx], -1]

        tic()
        rhipp_pme_aug_list[[time_idx]] <- pme(temp_rhipp, d = 2)
        temp_pme_aug_time <- toc()
        rhipp_pme_aug_time_list[[time_idx]] <- temp_pme_aug_time$toc -
          temp_pme_aug_time$tic
        rhipp_pme_aug_reconstruction_list[[
          time_idx
        ]] <- calculate_pme_reconstructions(
          rhipp_pme_aug_list[[time_idx]],
          temp_rhipp
        )
        rhipp_pme_aug_reconstruction_list[[time_idx]] <- cbind(
          time_values[time_idx],
          rhipp_pme_aug_reconstruction_list[[time_idx]]
        )

        tic()
        rhipp_pme_part[[1]][[time_idx]] <- pme(temp_rhipp_pt1, d = 2)
        temp_pme_time_pt1 <- toc()
        rhipp_pme_part_times[[1]][[time_idx]] <- temp_pme_time_pt1$toc -
          temp_pme_time_pt1$tic
        rhipp_pme_part_reconstructions[[1]][[
          time_idx
        ]] <- calculate_pme_reconstructions(
          rhipp_pme_part[[1]][[time_idx]],
          temp_rhipp_pt1[temp_rhipp_pt1[, 3] > 0, ]
        )
        tic()
        rhipp_pme_part[[2]][[time_idx]] <- pme(temp_rhipp_pt2, d = 2)
        temp_pme_time_pt2 <- toc()
        rhipp_pme_part_times[[2]][[time_idx]] <- temp_pme_time_pt2$toc -
          temp_pme_time_pt2$tic
        rhipp_pme_part_reconstructions[[2]][[
          time_idx
        ]] <- calculate_pme_reconstructions(
          rhipp_pme_part[[2]][[time_idx]],
          temp_rhipp_pt2[temp_rhipp_pt2[, 3] <= 0, ]
        )

        rhipp_pme_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          rhipp_pme_part_reconstructions[[1]][[time_idx]]
        )
        rhipp_pme_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          rhipp_pme_part_reconstructions[[2]][[time_idx]]
        )

        tic()
        principal_surface_part1 <- prinSurf(temp_rhipp_pt1)
        temp_pc_time_pt1 <- toc()
        rhipp_pc_part_times[[1]][[time_idx]] <- temp_pc_time_pt1$toc -
          temp_pc_time_pt1$tic
        surface_mse_part1 <- map(
          seq_along(principal_surface_part1),
          ~ principal_surface_part1[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part1 <- which.min(surface_mse_part1)

        rhipp_pc_part[[1]][[time_idx]] <- principal_surface_part1[[
          opt_surface_part1 + 2
        ]]
        rhipp_pc_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part1[[opt_surface_part1 + 2]]$PS
        )
        rhipp_pc_part_reconstructions[[1]][[
          time_idx
        ]] <- rhipp_pc_part_reconstructions[[1]][[time_idx]][
          temp_rhipp_pt1[, 3] > 0,
        ]

        tic()
        principal_surface_part2 <- prinSurf(temp_rhipp_pt2)
        temp_pc_time_pt2 <- toc()
        rhipp_pc_part_times[[2]][[time_idx]] <- temp_pc_time_pt2$toc -
          temp_pc_time_pt2$tic
        surface_mse_part2 <- map(
          seq_along(principal_surface_part2),
          ~ principal_surface_part2[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part2 <- which.min(surface_mse_part2)

        rhipp_pc_part[[2]][[time_idx]] <- principal_surface_part2[[
          opt_surface_part2 + 2
        ]]
        rhipp_pc_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part2[[opt_surface_part2 + 2]]$PS
        )
        rhipp_pc_part_reconstructions[[2]][[
          time_idx
        ]] <- rhipp_pc_part_reconstructions[[2]][[time_idx]][
          temp_rhipp_pt2[, 3] <= 0,
        ]
      }

      rhipp_lpme_part_reconstructions <- reduce(
        rhipp_lpme_part_reconstructions,
        rbind
      )

      rhipp_lpme_part_times <- reduce(
        rhipp_lpme_part_times,
        sum
      )

      rhipp_pme_aug_reconstructions <- reduce(
        rhipp_pme_aug_reconstruction_list,
        rbind
      )

      rhipp_pme_aug_time_list <- reduce(
        rhipp_pme_aug_time_list,
        sum
      )

      rhipp_pme_part_reconstructions <- map(
        seq_along(rhipp_pme_part_reconstructions),
        ~ reduce(rhipp_pme_part_reconstructions[[.x]], rbind)
      )
      rhipp_pme_part_reconstructions <- reduce(
        rhipp_pme_part_reconstructions,
        rbind
      )

      rhipp_pme_part_times <- map(
        seq_along(rhipp_pme_part_times),
        ~ reduce(rhipp_pme_part_times[[.x]], sum)
      )
      rhipp_pme_part_times <- reduce(
        rhipp_pme_part_times,
        sum
      )

      rhipp_pc_part_reconstructions <- map(
        seq_along(rhipp_pc_part_reconstructions),
        ~ reduce(rhipp_pc_part_reconstructions[[.x]], rbind)
      )
      rhipp_pc_part_reconstructions <- reduce(
        rhipp_pc_part_reconstructions,
        rbind
      )

      rhipp_pc_part_times <- map(
        seq_along(rhipp_pc_part_times),
        ~ reduce(rhipp_pc_part_times[[.x]], sum)
      )
      rhipp_pc_part_times <- reduce(
        rhipp_pc_part_times,
        sum
      )

      rhipp_lpme_part_volumes <- estimate_volume_interior_lpme(
        rhipp_lpme_part,
        list(rhipp_pt1, rhipp_pt2),
        time_values,
        n_points = 10000,
        data_max = patno_rhipp_centers,
        limit_scaler = 0.05
      )
      rhipp_lpme_part_volumes <- rhipp_lpme_part_volumes$volumes

      rhipp_pme_part_volumes <- estimate_volume_interior_pme(
        rhipp_pme_part,
        list(rhipp_pt1, rhipp_pt2),
        time_values,
        n_points = 10000,
        data_max = patno_rhipp_centers,
        limit_scaler = 0.05
      )
      rhipp_pme_part_volumes <- rhipp_pme_part_volumes$volumes

      rhipp_lpme_aug_volumes <- vector()
      rhipp_pme_aug_volumes <- vector()
      rhipp_pc_part_volumes <- vector()

      rhipp_lpme_aug_interior <- list()
      rhipp_pme_aug_interior <- list()
      rhipp_pc_part_interior <- list()

      rhipp_lpme_part_volumes_mesh <- vector()
      rhipp_pme_part_volumes_mesh <- vector()

      for (time_idx in seq_along(time_values)) {
        x_scale <- patno_rhipp_centers$max_x[time_idx]
        y_scale <- patno_rhipp_centers$max_y[time_idx]
        z_scale <- patno_rhipp_centers$max_z[time_idx]
        temp_lpme_reconstructions <- rhipp_lpme_aug_reconstructions[
          rhipp_lpme_aug_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pme_reconstructions <- rhipp_pme_aug_reconstructions[
          rhipp_pme_aug_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pc_reconstructions <- rhipp_pc_part_reconstructions[
          rhipp_pc_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]

        temp_lpme_part_reconstructions <- rhipp_lpme_part_reconstructions[
          rhipp_lpme_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pme_part_reconstructions <- rhipp_pme_part_reconstructions[
          rhipp_pme_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]

        temp_lpme_ashape <- ashape3d(
          temp_lpme_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        rhipp_lpme_aug_volumes[time_idx] <- volume_ashape3d(
          temp_lpme_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pme_ashape <- ashape3d(
          temp_pme_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        rhipp_pme_aug_volumes[time_idx] <- volume_ashape3d(
          temp_pme_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pc_ashape <- ashape3d(
          temp_pc_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        rhipp_pc_part_volumes[time_idx] <- volume_ashape3d(
          temp_pc_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_lpme_part_ashape <- ashape3d(
          temp_lpme_part_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        rhipp_lpme_part_volumes_mesh[time_idx] <- volume_ashape3d(
          temp_lpme_part_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pme_part_ashape <- ashape3d(
          temp_pme_part_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        rhipp_pme_part_volumes_mesh[time_idx] <- volume_ashape3d(
          temp_pme_part_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)
      }

      rhipp_lpme_aug_out <- list(
        lpme = rhipp_lpme_aug,
        reconstructions = rhipp_lpme_aug_reconstructions,
        volumes = rhipp_lpme_aug_volumes,
        fit_time = rhipp_lpme_aug_time
      )
      rhipp_pme_aug_out <- list(
        pme = rhipp_pme_aug_list,
        reconstructions = rhipp_pme_aug_reconstructions,
        volumes = rhipp_pme_aug_volumes,
        fit_time = rhipp_pme_aug_time_list
      )

      rhipp_lpme_part_out <- list(
        lpme = rhipp_lpme_part,
        reconstructions = rhipp_lpme_part_reconstructions,
        volumes = rhipp_lpme_part_volumes,
        volumes_mesh = rhipp_lpme_part_volumes_mesh,
        fit_time = rhipp_lpme_part_times
      )

      rhipp_pme_part_out <- list(
        pme = rhipp_pme_part,
        reconstructions = rhipp_pme_part_reconstructions,
        volumes = rhipp_pme_part_volumes,
        volumes_mesh = rhipp_pme_part_volumes_mesh,
        fit_time = rhipp_pme_part_times
      )

      rhipp_pc_part_out <- list(
        pc = rhipp_pc_part,
        reconstructions = rhipp_pc_part_reconstructions,
        volumes = rhipp_pc_part_volumes,
        fit_time = rhipp_pc_part_times
      )

      models <- list(
        lpme_aug = rhipp_lpme_aug_out,
        pme_aug = rhipp_pme_aug_out,
        lpme_part = rhipp_lpme_part_out,
        pme_part = rhipp_pme_part_out,
        pc_part = rhipp_pc_part_out
      )

      if (!dir.exists(here(paste0("output/adni/", patno)))) {
        dir.create(here(paste0("output/adni/", patno)), recursive = TRUE)
      }

      saveRDS(
        models,
        here(paste0("output/adni/", patno, "/rhipp_results.RDS"))
      )

      models
    },
    seed = TRUE,
    patnos = patnos,
    patno_idx = patno_idx,
    rhipp_surface = rhipp_surface,
    rhipp_centers = rhipp_centers
  )
}

# rhipp_results_out <- purrr::map(rhipp_results, value)
rhipp_results_out <- purrr::map(rhipp_results, ~ .x[])



rhipp_df_list <- list()

for (patno_idx in seq_along(patnos)) {
  patno <- patnos[patno_idx]
  rhipp <- rhipp_results_out[[patno_idx]]

  rhipp_lpme_aug <- rhipp$lpme_aug
  rhipp_pme_aug <- rhipp$pme_aug
  rhipp_lpme_part <- rhipp$lpme_part
  rhipp_pme_part <- rhipp$pme_part
  rhipp_pc_part <- rhipp$pc_part

  rhipp_time_values <- unique(rhipp_lpme_aug$reconstructions[, 1])
  rhipp_lpme_aug_vols <- rhipp_lpme_aug$volumes$volumes
  rhipp_lpme_aug_time <- rhipp_lpme_aug$fit_time

  rhipp_pme_aug_vols <- rhipp_pme_aug$volumes$volumes
  rhipp_pme_aug_time <- rhipp_pme_aug$fit_time

  rhipp_lpme_part_vols <- rhipp_lpme_part$volumes$volumes
  rhipp_lpme_part_vols_mesh <- rhipp_lpme_part$volumes_mesh
  rhipp_lpme_part_time <- rhipp_lpme_part$fit_time

  rhipp_pme_part_vols <- rhipp_pme_part$volumes$volumes
  rhipp_pme_part_vols_mesh <- rhipp_pme_part$volumes_mesh
  rhipp_pme_part_time <- rhipp_pme_part$fit_time

  rhipp_pc_part_vols <- rhipp_pc_part$volumes$volumes
  rhipp_pc_part_time <- rhipp_pc_part$fit_time

  rhipp_df_list[[patno_idx]] <- tibble(
    patno = patno,
    rhipp_time_values = rhipp_time_values,
    rhipp_lpme_aug_vols = rhipp_lpme_aug_vols,
    rhipp_pme_aug_vols = rhipp_pme_aug_vols,
    rhipp_lpme_part_vols = rhipp_lpme_part_vols,
    rhipp_pme_part_vols = rhipp_pme_part_vols,
    rhipp_pc_part_vols = rhipp_pc_part_vols,
    rhipp_lpme_part_vols_mesh = rhipp_lpme_part_vols_mesh,
    rhipp_pme_part_vols_mesh = rhipp_pme_part_vols_mesh,
    rhipp_lpme_aug_time = rhipp_lpme_aug_time,
    rhipp_pme_aug_time = rhipp_pme_aug_time,
    rhipp_lpme_part_time = rhipp_lpme_part_time,
    rhipp_pme_part_time = rhipp_pme_part_time,
    rhipp_pc_part_time = rhipp_pc_part_time
  )
}

rhipp_df <- bind_rows(rhipp_df_list)

lthal_results <- list()

for (patno_idx in seq_along(patnos)) {
  lthal_results[[patno_idx]] <- mirai(
    {
      patno <- patnos[patno_idx]
      patno_lthal <- lthal_surface |>
        filter(subid == patno) |>
        select(time_from_bl, x, y, z, theta, phi) |>
        mutate(
          partition = z > 0,
          partition1 = z > -0.15,
          partition2 = z < 0.15
        )

      patno_lthal_centers <- lthal_centers |>
        filter(subid == patno) |>
        group_by(date) |>
        summarize(
          max_x = max(max_x),
          max_y = max(max_y),
          max_z = max(max_z)
        )

      lthal_aug <- patno_lthal |>
        select(-partition, -partition1, -partition2) |>
        as.matrix()

      lthal_pt1 <- patno_lthal |>
        filter(partition1 == TRUE) |>
        select(time_from_bl, x, y, z) |>
        as.matrix()
      lthal_pt2 <- patno_lthal |>
        filter(partition2 == TRUE) |>
        select(time_from_bl, x, y, z) |>
        as.matrix()

      time_values <- unique(patno_lthal$time_from_bl)

      lthal_pme_aug_list <- list()
      lthal_pme_aug_time_list <- list()
      lthal_pme_aug_reconstruction_list <- list()

      lthal_lpme_part <- replicate(2, NULL, simplify = FALSE)
      lthal_pme_part <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      lthal_pc_part <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      lthal_lpme_part_reconstructions <- replicate(2, NULL, simplify = FALSE)
      lthal_pme_part_reconstructions <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      lthal_pc_part_reconstructions <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      lthal_lpme_part_times <- replicate(2, NULL, simplify = FALSE)
      lthal_pme_part_times <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      lthal_pc_part_times <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      tic()
      lthal_lpme_aug <- lpme(
        lthal_aug,
        d = 2,
        verbose = TRUE,
        print_plots = FALSE
      )
      lthal_lpme_aug_time <- toc()
      lthal_lpme_aug_time <- lthal_lpme_aug_time$toc - lthal_lpme_aug_time$tic

      lthal_lpme_aug_reconstructions <- calculate_lpme_reconstructions(
        lthal_lpme_aug,
        lthal_aug
      )

      tic()
      lthal_lpme_part[[1]] <- lpme(
        lthal_pt1,
        d = 2,
        verbose = FALSE,
        print_plots = FALSE
      )
      lthal_lpme_time_pt1 <- toc()
      lthal_lpme_part_times[[1]] <- lthal_lpme_time_pt1$toc -
        lthal_lpme_time_pt1$tic
      lthal_lpme_part_reconstructions[[1]] <- calculate_lpme_reconstructions(
        lthal_lpme_part[[1]],
        lthal_pt1[lthal_pt1[, 4] > 0, ]
      )
      tic()
      lthal_lpme_part[[2]] <- lpme(
        lthal_pt2,
        d = 2,
        verbose = FALSE,
        print_plots = FALSE
      )
      lthal_lpme_time_pt2 <- toc()
      lthal_lpme_part_times[[2]] <- lthal_lpme_time_pt2$toc -
        lthal_lpme_time_pt2$tic
      lthal_lpme_part_reconstructions[[2]] <- calculate_lpme_reconstructions(
        lthal_lpme_part[[2]],
        lthal_pt2[lthal_pt2[, 4] <= 0, ]
      )

      for (time_idx in seq_along(time_values)) {
        temp_lthal <- lthal_aug[lthal_aug[, 1] == time_values[time_idx], -1]
        temp_lthal_pt1 <- lthal_pt1[lthal_pt1[, 1] == time_values[time_idx], -1]
        temp_lthal_pt2 <- lthal_pt2[lthal_pt2[, 1] == time_values[time_idx], -1]

        tic()
        lthal_pme_aug_list[[time_idx]] <- pme(temp_lthal, d = 2)
        temp_pme_aug_time <- toc()
        lthal_pme_aug_time_list[[time_idx]] <- temp_pme_aug_time$toc -
          temp_pme_aug_time$tic
        lthal_pme_aug_reconstruction_list[[
          time_idx
        ]] <- calculate_pme_reconstructions(
          lthal_pme_aug_list[[time_idx]],
          temp_lthal
        )
        lthal_pme_aug_reconstruction_list[[time_idx]] <- cbind(
          time_values[time_idx],
          lthal_pme_aug_reconstruction_list[[time_idx]]
        )

        tic()
        lthal_pme_part[[1]][[time_idx]] <- pme(temp_lthal_pt1, d = 2)
        temp_pme_time_pt1 <- toc()
        lthal_pme_part_times[[1]][[time_idx]] <- temp_pme_time_pt1$toc -
          temp_pme_time_pt1$tic
        lthal_pme_part_reconstructions[[1]][[
          time_idx
        ]] <- calculate_pme_reconstructions(
          lthal_pme_part[[1]][[time_idx]],
          temp_lthal_pt1[temp_lthal_pt1[, 3] > 0, ]
        )
        tic()
        lthal_pme_part[[2]][[time_idx]] <- pme(temp_lthal_pt2, d = 2)
        temp_pme_time_pt2 <- toc()
        lthal_pme_part_times[[2]][[time_idx]] <- temp_pme_time_pt2$toc -
          temp_pme_time_pt2$tic
        lthal_pme_part_reconstructions[[2]][[
          time_idx
        ]] <- calculate_pme_reconstructions(
          lthal_pme_part[[2]][[time_idx]],
          temp_lthal_pt2[temp_lthal_pt2[, 3] <= 0, ]
        )

        lthal_pme_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          lthal_pme_part_reconstructions[[1]][[time_idx]]
        )
        lthal_pme_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          lthal_pme_part_reconstructions[[2]][[time_idx]]
        )

        tic()
        principal_surface_part1 <- prinSurf(temp_lthal_pt1)
        temp_pc_time_pt1 <- toc()
        lthal_pc_part_times[[1]][[time_idx]] <- temp_pc_time_pt1$toc -
          temp_pc_time_pt1$tic
        surface_mse_part1 <- map(
          seq_along(principal_surface_part1),
          ~ principal_surface_part1[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part1 <- which.min(surface_mse_part1)

        lthal_pc_part[[1]][[time_idx]] <- principal_surface_part1[[
          opt_surface_part1 + 2
        ]]
        lthal_pc_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part1[[opt_surface_part1 + 2]]$PS
        )
        lthal_pc_part_reconstructions[[1]][[
          time_idx
        ]] <- lthal_pc_part_reconstructions[[1]][[time_idx]][
          temp_lthal_pt1[, 3] > 0,
        ]

        tic()
        principal_surface_part2 <- prinSurf(temp_lthal_pt2)
        temp_pc_time_pt2 <- toc()
        lthal_pc_part_times[[2]][[time_idx]] <- temp_pc_time_pt2$toc -
          temp_pc_time_pt2$tic
        surface_mse_part2 <- map(
          seq_along(principal_surface_part2),
          ~ principal_surface_part2[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part2 <- which.min(surface_mse_part2)

        lthal_pc_part[[2]][[time_idx]] <- principal_surface_part2[[
          opt_surface_part2 + 2
        ]]
        lthal_pc_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part2[[opt_surface_part2 + 2]]$PS
        )
        lthal_pc_part_reconstructions[[2]][[
          time_idx
        ]] <- lthal_pc_part_reconstructions[[2]][[time_idx]][
          temp_lthal_pt2[, 3] <= 0,
        ]
      }

      lthal_lpme_part_reconstructions <- reduce(
        lthal_lpme_part_reconstructions,
        rbind
      )

      lthal_lpme_part_times <- reduce(
        lthal_lpme_part_times,
        sum
      )

      lthal_pme_aug_reconstructions <- reduce(
        lthal_pme_aug_reconstruction_list,
        rbind
      )

      lthal_pme_aug_time_list <- reduce(
        lthal_pme_aug_time_list,
        sum
      )

      lthal_pme_part_reconstructions <- map(
        seq_along(lthal_pme_part_reconstructions),
        ~ reduce(lthal_pme_part_reconstructions[[.x]], rbind)
      )
      lthal_pme_part_reconstructions <- reduce(
        lthal_pme_part_reconstructions,
        rbind
      )

      lthal_pme_part_times <- map(
        seq_along(lthal_pme_part_times),
        ~ reduce(lthal_pme_part_times[[.x]], sum)
      )
      lthal_pme_part_times <- reduce(
        lthal_pme_part_times,
        sum
      )

      lthal_pc_part_reconstructions <- map(
        seq_along(lthal_pc_part_reconstructions),
        ~ reduce(lthal_pc_part_reconstructions[[.x]], rbind)
      )
      lthal_pc_part_reconstructions <- reduce(
        lthal_pc_part_reconstructions,
        rbind
      )

      lthal_pc_part_times <- map(
        seq_along(lthal_pc_part_times),
        ~ reduce(lthal_pc_part_times[[.x]], sum)
      )
      lthal_pc_part_times <- reduce(
        lthal_pc_part_times,
        sum
      )

      lthal_lpme_part_volumes <- estimate_volume_interior_lpme(
        lthal_lpme_part,
        list(lthal_pt1, lthal_pt2),
        time_values,
        n_points = 10000,
        data_max = patno_lthal_centers,
        limit_scaler = 0.05
      )
      lthal_lpme_part_volumes <- lthal_lpme_part_volumes$volumes

      lthal_pme_part_volumes <- estimate_volume_interior_pme(
        lthal_pme_part,
        list(lthal_pt1, lthal_pt2),
        time_values,
        n_points = 10000,
        data_max = patno_lthal_centers,
        limit_scaler = 0.05
      )
      lthal_pme_part_volumes <- lthal_pme_part_volumes$volumes

      lthal_lpme_aug_volumes <- vector()
      lthal_pme_aug_volumes <- vector()
      lthal_pc_part_volumes <- vector()

      lthal_lpme_aug_interior <- list()
      lthal_pme_aug_interior <- list()
      lthal_pc_part_interior <- list()

      lthal_lpme_part_volumes_mesh <- vector()
      lthal_pme_part_volumes_mesh <- vector()

      for (time_idx in seq_along(time_values)) {
        x_scale <- patno_lthal_centers$max_x[time_idx]
        y_scale <- patno_lthal_centers$max_y[time_idx]
        z_scale <- patno_lthal_centers$max_z[time_idx]
        temp_lpme_reconstructions <- lthal_lpme_aug_reconstructions[
          lthal_lpme_aug_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pme_reconstructions <- lthal_pme_aug_reconstructions[
          lthal_pme_aug_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pc_reconstructions <- lthal_pc_part_reconstructions[
          lthal_pc_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]

        temp_lpme_part_reconstructions <- lthal_lpme_part_reconstructions[
          lthal_lpme_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pme_part_reconstructions <- lthal_pme_part_reconstructions[
          lthal_pme_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]

        temp_lpme_ashape <- ashape3d(
          temp_lpme_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        lthal_lpme_aug_volumes[time_idx] <- volume_ashape3d(
          temp_lpme_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pme_ashape <- ashape3d(
          temp_pme_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        lthal_pme_aug_volumes[time_idx] <- volume_ashape3d(
          temp_pme_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pc_ashape <- ashape3d(
          temp_pc_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        lthal_pc_part_volumes[time_idx] <- volume_ashape3d(
          temp_pc_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_lpme_part_ashape <- ashape3d(
          temp_lpme_part_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        lthal_lpme_part_volumes_mesh[time_idx] <- volume_ashape3d(
          temp_lpme_part_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pme_part_ashape <- ashape3d(
          temp_pme_part_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        lthal_pme_part_volumes_mesh[time_idx] <- volume_ashape3d(
          temp_pme_part_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)
      }

      lthal_lpme_aug_out <- list(
        lpme = lthal_lpme_aug,
        reconstructions = lthal_lpme_aug_reconstructions,
        volumes = lthal_lpme_aug_volumes,
        fit_time = lthal_lpme_aug_time
      )
      lthal_pme_aug_out <- list(
        pme = lthal_pme_aug_list,
        reconstructions = lthal_pme_aug_reconstructions,
        volumes = lthal_pme_aug_volumes,
        fit_time = lthal_pme_aug_time_list
      )

      lthal_lpme_part_out <- list(
        lpme = lthal_lpme_part,
        reconstructions = lthal_lpme_part_reconstructions,
        volumes = lthal_lpme_part_volumes,
        volumes_mesh = lthal_lpme_part_volumes_mesh,
        fit_time = lthal_lpme_part_times
      )

      lthal_pme_part_out <- list(
        pme = lthal_pme_part,
        reconstructions = lthal_pme_part_reconstructions,
        volumes = lthal_pme_part_volumes,
        volumes_mesh = lthal_pme_part_volumes_mesh,
        fit_time = lthal_pme_part_times
      )

      lthal_pc_part_out <- list(
        pc = lthal_pc_part,
        reconstructions = lthal_pc_part_reconstructions,
        volumes = lthal_pc_part_volumes,
        fit_time = lthal_pc_part_times
      )

      models <- list(
        lpme_aug = lthal_lpme_aug_out,
        pme_aug = lthal_pme_aug_out,
        lpme_part = lthal_lpme_part_out,
        pme_part = lthal_pme_part_out,
        pc_part = lthal_pc_part_out
      )

      if (!dir.exists(here(paste0("output/adni/", patno)))) {
        dir.create(here(paste0("output/adni/", patno)), recursive = TRUE)
      }

      saveRDS(
        models,
        here(paste0("output/adni/", patno, "/lthal_results.RDS"))
      )

      models
    },
    seed = TRUE,
    patnos = patnos,
    patno_idx = patno_idx,
    lthal_surface = lthal_surface,
    lthal_centers = lthal_centers
  )
}

# lthal_results_out <- purrr::map(lthal_results, value)
lthal_results_out <- purrr::map(lthal_results, ~ .x[])



lthal_df_list <- list()

for (patno_idx in seq_along(patnos)) {
  patno <- patnos[patno_idx]
  lthal <- lthal_results_out[[patno_idx]]

  lthal_lpme_aug <- lthal$lpme_aug
  lthal_pme_aug <- lthal$pme_aug
  lthal_lpme_part <- lthal$lpme_part
  lthal_pme_part <- lthal$pme_part
  lthal_pc_part <- lthal$pc_part

  lthal_time_values <- unique(lthal_lpme_aug$reconstructions[, 1])
  lthal_lpme_aug_vols <- lthal_lpme_aug$volumes$volumes
  lthal_lpme_aug_time <- lthal_lpme_aug$fit_time

  lthal_pme_aug_vols <- lthal_pme_aug$volumes$volumes
  lthal_pme_aug_time <- lthal_pme_aug$fit_time

  lthal_lpme_part_vols <- lthal_lpme_part$volumes$volumes
  lthal_lpme_part_vols_mesh <- lthal_lpme_part$volumes_mesh
  lthal_lpme_part_time <- lthal_lpme_part$fit_time

  lthal_pme_part_vols <- lthal_pme_part$volumes$volumes
  lthal_pme_part_vols_mesh <- lthal_pme_part$volumes_mesh
  lthal_pme_part_time <- lthal_pme_part$fit_time

  lthal_pc_part_vols <- lthal_pc_part$volumes$volumes
  lthal_pc_part_time <- lthal_pc_part$fit_time

  lthal_df_list[[patno_idx]] <- tibble(
    patno = patno,
    lthal_time_values = lthal_time_values,
    lthal_lpme_aug_vols = lthal_lpme_aug_vols,
    lthal_pme_aug_vols = lthal_pme_aug_vols,
    lthal_lpme_part_vols = lthal_lpme_part_vols,
    lthal_pme_part_vols = lthal_pme_part_vols,
    lthal_pc_part_vols = lthal_pc_part_vols,
    lthal_lpme_part_vols_mesh = lthal_lpme_part_vols_mesh,
    lthal_pme_part_vols_mesh = lthal_pme_part_vols_mesh,
    lthal_lpme_aug_time = lthal_lpme_aug_time,
    lthal_pme_aug_time = lthal_pme_aug_time,
    lthal_lpme_part_time = lthal_lpme_part_time,
    lthal_pme_part_time = lthal_pme_part_time,
    lthal_pc_part_time = lthal_pc_part_time
  )
}

lthal_df <- bind_rows(lthal_df_list)

rthal_results <- list()

for (patno_idx in seq_along(patnos)) {
  rthal_results[[patno_idx]] <- mirai(
    {
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
        limit_scaler = 0.05
      )
      rthal_lpme_part_volumes <- rthal_lpme_part_volumes$volumes

      rthal_pme_part_volumes <- estimate_volume_interior_pme(
        rthal_pme_part,
        list(rthal_pt1, rthal_pt2),
        time_values,
        n_points = 10000,
        data_max = patno_rthal_centers,
        limit_scaler = 0.05
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

      models <- list(
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
        models,
        here(paste0("output/adni/", patno, "/rthal_results.RDS"))
      )

      models
    },
    seed = TRUE,
    patnos = patnos,
    patno_idx = patno_idx,
    rthal_surface = rthal_surface,
    rthal_centers = rthal_centers
  )
}

# rthal_results_out <- purrr::map(rthal_results, value)
rthal_results_out <- purrr::map(rthal_results, ~ .x[])



rthal_df_list <- list()

for (patno_idx in seq_along(patnos)) {
  patno <- patnos[patno_idx]
  rthal <- rthal_results_out[[patno_idx]]

  rthal_lpme_aug <- rthal$lpme_aug
  rthal_pme_aug <- rthal$pme_aug
  rthal_lpme_part <- rthal$lpme_part
  rthal_pme_part <- rthal$pme_part
  rthal_pc_part <- rthal$pc_part

  rthal_time_values <- unique(rthal_lpme_aug$reconstructions[, 1])
  rthal_lpme_aug_vols <- rthal_lpme_aug$volumes$volumes
  rthal_lpme_aug_time <- rthal_lpme_aug$fit_time

  rthal_pme_aug_vols <- rthal_pme_aug$volumes$volumes
  rthal_pme_aug_time <- rthal_pme_aug$fit_time

  rthal_lpme_part_vols <- rthal_lpme_part$volumes$volumes
  rthal_lpme_part_vols_mesh <- rthal_lpme_part$volumes_mesh
  rthal_lpme_part_time <- rthal_lpme_part$fit_time

  rthal_pme_part_vols <- rthal_pme_part$volumes$volumes
  rthal_pme_part_vols_mesh <- rthal_pme_part$volumes_mesh
  rthal_pme_part_time <- rthal_pme_part$fit_time

  rthal_pc_part_vols <- rthal_pc_part$volumes$volumes
  rthal_pc_part_time <- rthal_pc_part$fit_time

  rthal_df_list[[patno_idx]] <- tibble(
    patno = patno,
    rthal_time_values = rthal_time_values,
    rthal_lpme_aug_vols = rthal_lpme_aug_vols,
    rthal_pme_aug_vols = rthal_pme_aug_vols,
    rthal_lpme_part_vols = rthal_lpme_part_vols,
    rthal_pme_part_vols = rthal_pme_part_vols,
    rthal_pc_part_vols = rthal_pc_part_vols,
    rthal_lpme_part_vols_mesh = rthal_lpme_part_vols_mesh,
    rthal_pme_part_vols_mesh = rthal_pme_part_vols_mesh,
    rthal_lpme_aug_time = rthal_lpme_aug_time,
    rthal_pme_aug_time = rthal_pme_aug_time,
    rthal_lpme_part_time = rthal_lpme_part_time,
    rthal_pme_part_time = rthal_pme_part_time,
    rthal_pc_part_time = rthal_pc_part_time
  )
}

rthal_df <- bind_rows(rthal_df_list)

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
