library(here)
library(lubridate)
# library(mirai)
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
```


```{r}
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
    max_x = max(abs(x)),
    max_y = max(abs(y)),
    max_z = max(abs(z))
  ) |>
  ungroup()

rhipp_centers <- rhipp_surface |>
  group_by(subid, date, scan_id) |>
  summarise(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x)),
    max_y = max(abs(y)),
    max_z = max(abs(z))
  ) |>
  ungroup()

lthal_centers <- lthal_surface |>
  group_by(subid, date, scan_id) |>
  summarise(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x)),
    max_y = max(abs(y)),
    max_z = max(abs(z))
  ) |>
  ungroup()

rthal_centers <- rthal_surface |>
  group_by(subid, date, scan_id) |>
  summarise(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x)),
    max_y = max(abs(y)),
    max_z = max(abs(z))
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
patnos <- patnos[1:5]
cores <- parallel::detectCores()
plan(multisession, workers = cores)

lhipp_results <- list()

for (patno_idx in seq_along(patnos)) {
  lhipp_results[[patno_idx]] <- future(
    {
      patno <- patnos[patno_idx]
      patno_lhipp <- lhipp_surface |>
        filter(subid == patno) |>
        select(time_from_bl, x, y, z, theta, phi) |>
        mutate(
          partition = z > 0,
          partition1 = z > -0.025,
          partition2 = z < 0.025
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
        verbose = FALSE, 
        print_plots = FALSE
      )
      lhipp_lpme_aug_time <- toc()

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
      lhipp_lpme_part_times[[1]] <- toc()
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
      lhipp_lpme_part_times[[2]] <- toc()
      lhipp_lpme_part_reconstructions[[2]] <- calculate_lpme_reconstructions(
        lhipp_lpme_part[[2]],
        lhipp_pt2[lhipp_pt2[, 4] <= 0, ]
      )

      for (time_idx in seq_along(time_values)) {
        temp_lhipp <- lhipp_aug[lhipp_aug[, 1 ] == time_values[time_idx], -1]
        temp_lhipp_pt1 <- lhipp_pt1[lhipp_pt1[, 1] == time_values[time_idx], -1]
        temp_lhipp_pt2 <- lhipp_pt2[lhipp_pt2[, 1] == time_values[time_idx], -1]

        tic()
        lhipp_pme_aug_list[[time_idx]] <- pme(temp_lhipp, d = 2)
        lhipp_pme_aug_time_list[[time_idx]] <- toc()
        lhipp_pme_aug_reconstruction_list[[time_idx]] <- calculate_pme_reconstructions(
          lhipp_pme_aug_list[[time_idx]],
          temp_lhipp
        )
        lhipp_pme_aug_reconstruction_list[[time_idx]] <- cbind(
          time_values[time_idx],
          lhipp_pme_aug_reconstruction_list[[time_idx]]
        )


        tic()
        lhipp_pme_part[[1]][[time_idx]] <- pme(temp_lhipp_pt1, d = 2)
        lhipp_pme_part_times[[1]][[time_idx]] <- toc()
        lhipp_pme_part_reconstructions[[1]][[time_idx]] <- calculate_pme_reconstructions(
          lhipp_pme_part[[1]][[time_idx]],
          temp_lhipp_pt1[temp_lhipp_pt1[, 3] > 0, ]
        )
        tic()
        lhipp_pme_part[[2]][[time_idx]] <- pme(temp_lhipp_pt2, d = 2)
        lhipp_pme_part_times[[2]][[time_idx]] <- toc()
        lhipp_pme_part_reconstructions[[2]][[time_idx]] <- calculate_pme_reconstructions(
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
        lhipp_pc_part_times[[1]][[time_idx]] <- toc()
        surface_mse_part1 <- map(
          seq_along(principal_surface_part1),
          ~ principal_surface_part1[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part1 <- which.min(surface_mse_part1)

        lhipp_pc_part[[1]][[time_idx]] <- principal_surface_part1[[opt_surface_part1 + 2]]
        lhipp_pc_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part1[[opt_surface_part1 + 2]]$PS
        )
        lhipp_pc_part_reconstructions[[1]][[time_idx]] <- lhipp_pc_part_reconstructions[[1]][[time_idx]][temp_lhipp_pt1[, 3] > 0, ]

        tic()
        principal_surface_part2 <- prinSurf(temp_lhipp_pt2)
        lhipp_pc_part_times[[2]][[time_idx]] <- toc()
        surface_mse_part2 <- map(
          seq_along(principal_surface_part2),
          ~ principal_surface_part2[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part2 <- which.min(surface_mse_part2)

        lhipp_pc_part[[2]][[time_idx]] <- principal_surface_part2[[opt_surface_part2 + 2]]
        lhipp_pc_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part2[[opt_surface_part2 + 2]]$PS
        )
        lhipp_pc_part_reconstructions[[2]][[time_idx]] <- lhipp_pc_part_reconstructions[[2]][[time_idx]][temp_lhipp_pt2[, 3] <= 0, ]
      }

      lhipp_lpme_part_reconstructions <- reduce(
        lhipp_lpme_part_reconstructions,
        rbind
      )

      lhipp_pme_aug_reconstructions  <- reduce(
        lhipp_pme_aug_reconstruction_list,
        rbind
      )

      lhipp_pme_part_reconstructions <- map(
        seq_along(lhipp_pme_part_reconstructions),
        ~ reduce(lhipp_pme_part_reconstructions[[.x]], rbind)
      )
      lhipp_pme_part_reconstructions <- reduce(
        lhipp_pme_part_reconstructions, 
        rbind
      )

      lhipp_pc_part_reconstructions <- map(
        seq_along(lhipp_pc_part_reconstructions),
        ~ reduce(lhipp_pc_part_reconstructions[[.x]], rbind)
      )
      lhipp_pc_part_reconstructions <- reduce(
        lhipp_pc_part_reconstructions, 
        rbind
      )

      lhipp_lpme_aug_out <- list(
        lpme = lhipp_lpme_aug,
        reconstructions = lhipp_lpme_aug_reconstructions,
        fit_time = lhipp_lpme_aug_time
      )
      lhipp_pme_aug_out <- list(
        pme = lhipp_pme_aug_list,
        reconstructions = lhipp_pme_aug_reconstructions,
        fit_time = lhipp_pme_aug_time_list
      )

      lhipp_lpme_part_out <- list(
        lpme = lhipp_lpme_part,
        reconstructions = lhipp_lpme_part_reconstructions,
        fit_time = lhipp_lpme_part_times
      )

      lhipp_pme_part_out <- list(
        pme = lhipp_pme_part,
        reconstructions = lhipp_pme_part_reconstructions,
        fit_time = lhipp_pme_part_times
      )

      lhipp_pc_part_out <- list(
        pc = lhipp_pc_part,
        reconstructions = lhipp_pc_part_reconstructions,
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
    seed = TRUE
    # patno_idx = patno_idx,
    # patnos = patnos,
    # lhipp_surface = lhipp_surface
  )
}

lhipp_results <- map(lhipp_results, ~ value(.x))


plan(multisession, workers = cores)

rhipp_results <- list()

for (patno_idx in seq_along(patnos)) {
  rhipp_results[[patno_idx]] <- future(
    {
      patno <- patnos[patno_idx]
      patno_rhipp <- rhipp_surface |>
        filter(subid == patno) |>
        select(time_from_bl, x, y, z, theta, phi) |>
        mutate(
          partition = z > 0,
          partition1 = z > -0.025,
          partition2 = z < 0.025
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
        verbose = FALSE, 
        print_plots = FALSE
      )
      rhipp_lpme_aug_time <- toc()

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
      rhipp_lpme_part_times[[1]] <- toc()
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
      rhipp_lpme_part_times[[2]] <- toc()
      rhipp_lpme_part_reconstructions[[2]] <- calculate_lpme_reconstructions(
        rhipp_lpme_part[[2]],
        rhipp_pt2[rhipp_pt2[, 4] <= 0, ]
      )

      for (time_idx in seq_along(time_values)) {
        temp_rhipp <- rhipp_aug[rhipp_aug[, 1 ] == time_values[time_idx], -1]
        temp_rhipp_pt1 <- rhipp_pt1[rhipp_pt1[, 1] == time_values[time_idx], -1]
        temp_rhipp_pt2 <- rhipp_pt2[rhipp_pt2[, 1] == time_values[time_idx], -1]

        tic()
        rhipp_pme_aug_list[[time_idx]] <- pme(temp_rhipp, d = 2)
        rhipp_pme_aug_time_list[[time_idx]] <- toc()
        rhipp_pme_aug_reconstruction_list[[time_idx]] <- calculate_pme_reconstructions(
          rhipp_pme_aug_list[[time_idx]],
          temp_rhipp
        )
        rhipp_pme_aug_reconstruction_list[[time_idx]] <- cbind(
          time_values[time_idx],
          rhipp_pme_aug_reconstruction_list[[time_idx]]
        )


        tic()
        rhipp_pme_part[[1]][[time_idx]] <- pme(temp_rhipp_pt1, d = 2)
        rhipp_pme_part_times[[1]][[time_idx]] <- toc()
        rhipp_pme_part_reconstructions[[1]][[time_idx]] <- calculate_pme_reconstructions(
          rhipp_pme_part[[1]][[time_idx]],
          temp_rhipp_pt1[temp_rhipp_pt1[, 3] > 0, ]
        )
        tic()
        rhipp_pme_part[[2]][[time_idx]] <- pme(temp_rhipp_pt2, d = 2)
        rhipp_pme_part_times[[2]][[time_idx]] <- toc()
        rhipp_pme_part_reconstructions[[2]][[time_idx]] <- calculate_pme_reconstructions(
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
        rhipp_pc_part_times[[1]][[time_idx]] <- toc()
        surface_mse_part1 <- map(
          seq_along(principal_surface_part1),
          ~ principal_surface_part1[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part1 <- which.min(surface_mse_part1)

        rhipp_pc_part[[1]][[time_idx]] <- principal_surface_part1[[opt_surface_part1 + 2]]
        rhipp_pc_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part1[[opt_surface_part1 + 2]]$PS
        )
        rhipp_pc_part_reconstructions[[1]][[time_idx]] <- rhipp_pc_part_reconstructions[[1]][[time_idx]][temp_rhipp_pt1[, 3] > 0, ]

        tic()
        principal_surface_part2 <- prinSurf(temp_rhipp_pt2)
        rhipp_pc_part_times[[2]][[time_idx]] <- toc()
        surface_mse_part2 <- map(
          seq_along(principal_surface_part2),
          ~ principal_surface_part2[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part2 <- which.min(surface_mse_part2)

        rhipp_pc_part[[2]][[time_idx]] <- principal_surface_part2[[opt_surface_part2 + 2]]
        rhipp_pc_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part2[[opt_surface_part2 + 2]]$PS
        )
        rhipp_pc_part_reconstructions[[2]][[time_idx]] <- rhipp_pc_part_reconstructions[[2]][[time_idx]][temp_rhipp_pt2[, 3] <= 0, ]
      }

      rhipp_lpme_part_reconstructions <- reduce(
        rhipp_lpme_part_reconstructions,
        rbind
      )

      rhipp_pme_aug_reconstructions  <- reduce(
        rhipp_pme_aug_reconstruction_list,
        rbind
      )

      rhipp_pme_part_reconstructions <- map(
        seq_along(rhipp_pme_part_reconstructions),
        ~ reduce(rhipp_pme_part_reconstructions[[.x]], rbind)
      )
      rhipp_pme_part_reconstructions <- reduce(
        rhipp_pme_part_reconstructions, 
        rbind
      )

      rhipp_pc_part_reconstructions <- map(
        seq_along(rhipp_pc_part_reconstructions),
        ~ reduce(rhipp_pc_part_reconstructions[[.x]], rbind)
      )
      rhipp_pc_part_reconstructions <- reduce(
        rhipp_pc_part_reconstructions, 
        rbind
      )

      rhipp_lpme_aug_out <- list(
        lpme = rhipp_lpme_aug,
        reconstructions = rhipp_lpme_aug_reconstructions,
        fit_time = rhipp_lpme_aug_time
      )
      rhipp_pme_aug_out <- list(
        pme = rhipp_pme_aug_list,
        reconstructions = rhipp_pme_aug_reconstructions,
        fit_time = rhipp_pme_aug_time_list
      )

      rhipp_lpme_part_out <- list(
        lpme = rhipp_lpme_part,
        reconstructions = rhipp_lpme_part_reconstructions,
        fit_time = rhipp_lpme_part_times
      )

      rhipp_pme_part_out <- list(
        pme = rhipp_pme_part,
        reconstructions = rhipp_pme_part_reconstructions,
        fit_time = rhipp_pme_part_times
      )

      rhipp_pc_part_out <- list(
        pc = rhipp_pc_part,
        reconstructions = rhipp_pc_part_reconstructions,
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
    seed = TRUE
    # patno_idx = patno_idx,
    # patnos = patnos,
    # rhipp_surface = rhipp_surface
  )
}

rhipp_results <- map(rhipp_results, ~ value(.x))


plan(multisession, workers = cores)

lthal_results <- list()

for (patno_idx in seq_along(patnos)) {
  lthal_results[[patno_idx]] <- future(
    {
      patno <- patnos[patno_idx]
      patno_lthal <- lthal_surface |>
        filter(subid == patno) |>
        select(time_from_bl, x, y, z, theta, phi) |>
        mutate(
          partition = z > 0,
          partition1 = z > -0.025,
          partition2 = z < 0.025
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
        verbose = FALSE, 
        print_plots = FALSE
      )
      lthal_lpme_aug_time <- toc()

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
      lthal_lpme_part_times[[1]] <- toc()
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
      lthal_lpme_part_times[[2]] <- toc()
      lthal_lpme_part_reconstructions[[2]] <- calculate_lpme_reconstructions(
        lthal_lpme_part[[2]],
        lthal_pt2[lthal_pt2[, 4] <= 0, ]
      )

      for (time_idx in seq_along(time_values)) {
        temp_lthal <- lthal_aug[lthal_aug[, 1 ] == time_values[time_idx], -1]
        temp_lthal_pt1 <- lthal_pt1[lthal_pt1[, 1] == time_values[time_idx], -1]
        temp_lthal_pt2 <- lthal_pt2[lthal_pt2[, 1] == time_values[time_idx], -1]

        tic()
        lthal_pme_aug_list[[time_idx]] <- pme(temp_lthal, d = 2)
        lthal_pme_aug_time_list[[time_idx]] <- toc()
        lthal_pme_aug_reconstruction_list[[time_idx]] <- calculate_pme_reconstructions(
          lthal_pme_aug_list[[time_idx]],
          temp_lthal
        )
        lthal_pme_aug_reconstruction_list[[time_idx]] <- cbind(
          time_values[time_idx],
          lthal_pme_aug_reconstruction_list[[time_idx]]
        )


        tic()
        lthal_pme_part[[1]][[time_idx]] <- pme(temp_lthal_pt1, d = 2)
        lthal_pme_part_times[[1]][[time_idx]] <- toc()
        lthal_pme_part_reconstructions[[1]][[time_idx]] <- calculate_pme_reconstructions(
          lthal_pme_part[[1]][[time_idx]],
          temp_lthal_pt1[temp_lthal_pt1[, 3] > 0, ]
        )
        tic()
        lthal_pme_part[[2]][[time_idx]] <- pme(temp_lthal_pt2, d = 2)
        lthal_pme_part_times[[2]][[time_idx]] <- toc()
        lthal_pme_part_reconstructions[[2]][[time_idx]] <- calculate_pme_reconstructions(
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
        lthal_pc_part_times[[1]][[time_idx]] <- toc()
        surface_mse_part1 <- map(
          seq_along(principal_surface_part1),
          ~ principal_surface_part1[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part1 <- which.min(surface_mse_part1)

        lthal_pc_part[[1]][[time_idx]] <- principal_surface_part1[[opt_surface_part1 + 2]]
        lthal_pc_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part1[[opt_surface_part1 + 2]]$PS
        )
        lthal_pc_part_reconstructions[[1]][[time_idx]] <- lthal_pc_part_reconstructions[[1]][[time_idx]][temp_lthal_pt1[, 3] > 0, ]

        tic()
        principal_surface_part2 <- prinSurf(temp_lthal_pt2)
        lthal_pc_part_times[[2]][[time_idx]] <- toc()
        surface_mse_part2 <- map(
          seq_along(principal_surface_part2),
          ~ principal_surface_part2[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part2 <- which.min(surface_mse_part2)

        lthal_pc_part[[2]][[time_idx]] <- principal_surface_part2[[opt_surface_part2 + 2]]
        lthal_pc_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part2[[opt_surface_part2 + 2]]$PS
        )
        lthal_pc_part_reconstructions[[2]][[time_idx]] <- lthal_pc_part_reconstructions[[2]][[time_idx]][temp_lthal_pt2[, 3] <= 0, ]
      }

      lthal_lpme_part_reconstructions <- reduce(
        lthal_lpme_part_reconstructions,
        rbind
      )

      lthal_pme_aug_reconstructions  <- reduce(
        lthal_pme_aug_reconstruction_list,
        rbind
      )

      lthal_pme_part_reconstructions <- map(
        seq_along(lthal_pme_part_reconstructions),
        ~ reduce(lthal_pme_part_reconstructions[[.x]], rbind)
      )
      lthal_pme_part_reconstructions <- reduce(
        lthal_pme_part_reconstructions, 
        rbind
      )

      lthal_pc_part_reconstructions <- map(
        seq_along(lthal_pc_part_reconstructions),
        ~ reduce(lthal_pc_part_reconstructions[[.x]], rbind)
      )
      lthal_pc_part_reconstructions <- reduce(
        lthal_pc_part_reconstructions, 
        rbind
      )

      lthal_lpme_aug_out <- list(
        lpme = lthal_lpme_aug,
        reconstructions = lthal_lpme_aug_reconstructions,
        fit_time = lthal_lpme_aug_time
      )
      lthal_pme_aug_out <- list(
        pme = lthal_pme_aug_list,
        reconstructions = lthal_pme_aug_reconstructions,
        fit_time = lthal_pme_aug_time_list
      )

      lthal_lpme_part_out <- list(
        lpme = lthal_lpme_part,
        reconstructions = lthal_lpme_part_reconstructions,
        fit_time = lthal_lpme_part_times
      )

      lthal_pme_part_out <- list(
        pme = lthal_pme_part,
        reconstructions = lthal_pme_part_reconstructions,
        fit_time = lthal_pme_part_times
      )

      lthal_pc_part_out <- list(
        pc = lthal_pc_part,
        reconstructions = lthal_pc_part_reconstructions,
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
    seed = TRUE
    # patno_idx = patno_idx,
    # patnos = patnos,
    # lthal_surface = lthal_surface
  )
}

lthal_results <- map(lthal_results, ~ value(.x))


plan(multisession, workers = cores)

rthal_results <- list()

for (patno_idx in seq_along(patnos)) {
  rthal_results[[patno_idx]] <- future(
    {
      patno <- patnos[patno_idx]
      patno_rthal <- rthal_surface |>
        filter(subid == patno) |>
        select(time_from_bl, x, y, z, theta, phi) |>
        mutate(
          partition = z > 0,
          partition1 = z > -0.025,
          partition2 = z < 0.025
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
        verbose = FALSE, 
        print_plots = FALSE
      )
      rthal_lpme_aug_time <- toc()

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
      rthal_lpme_part_times[[1]] <- toc()
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
      rthal_lpme_part_times[[2]] <- toc()
      rthal_lpme_part_reconstructions[[2]] <- calculate_lpme_reconstructions(
        rthal_lpme_part[[2]],
        rthal_pt2[rthal_pt2[, 4] <= 0, ]
      )

      for (time_idx in seq_along(time_values)) {
        temp_rthal <- rthal_aug[rthal_aug[, 1 ] == time_values[time_idx], -1]
        temp_rthal_pt1 <- rthal_pt1[rthal_pt1[, 1] == time_values[time_idx], -1]
        temp_rthal_pt2 <- rthal_pt2[rthal_pt2[, 1] == time_values[time_idx], -1]

        tic()
        rthal_pme_aug_list[[time_idx]] <- pme(temp_rthal, d = 2)
        rthal_pme_aug_time_list[[time_idx]] <- toc()
        rthal_pme_aug_reconstruction_list[[time_idx]] <- calculate_pme_reconstructions(
          rthal_pme_aug_list[[time_idx]],
          temp_rthal
        )
        rthal_pme_aug_reconstruction_list[[time_idx]] <- cbind(
          time_values[time_idx],
          rthal_pme_aug_reconstruction_list[[time_idx]]
        )


        tic()
        rthal_pme_part[[1]][[time_idx]] <- pme(temp_rthal_pt1, d = 2)
        rthal_pme_part_times[[1]][[time_idx]] <- toc()
        rthal_pme_part_reconstructions[[1]][[time_idx]] <- calculate_pme_reconstructions(
          rthal_pme_part[[1]][[time_idx]],
          temp_rthal_pt1[temp_rthal_pt1[, 3] > 0, ]
        )
        tic()
        rthal_pme_part[[2]][[time_idx]] <- pme(temp_rthal_pt2, d = 2)
        rthal_pme_part_times[[2]][[time_idx]] <- toc()
        rthal_pme_part_reconstructions[[2]][[time_idx]] <- calculate_pme_reconstructions(
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
        rthal_pc_part_times[[1]][[time_idx]] <- toc()
        surface_mse_part1 <- map(
          seq_along(principal_surface_part1),
          ~ principal_surface_part1[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part1 <- which.min(surface_mse_part1)

        rthal_pc_part[[1]][[time_idx]] <- principal_surface_part1[[opt_surface_part1 + 2]]
        rthal_pc_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part1[[opt_surface_part1 + 2]]$PS
        )
        rthal_pc_part_reconstructions[[1]][[time_idx]] <- rthal_pc_part_reconstructions[[1]][[time_idx]][temp_rthal_pt1[, 3] > 0, ]

        tic()
        principal_surface_part2 <- prinSurf(temp_rthal_pt2)
        rthal_pc_part_times[[2]][[time_idx]] <- toc()
        surface_mse_part2 <- map(
          seq_along(principal_surface_part2),
          ~ principal_surface_part2[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part2 <- which.min(surface_mse_part2)

        rthal_pc_part[[2]][[time_idx]] <- principal_surface_part2[[opt_surface_part2 + 2]]
        rthal_pc_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part2[[opt_surface_part2 + 2]]$PS
        )
        rthal_pc_part_reconstructions[[2]][[time_idx]] <- rthal_pc_part_reconstructions[[2]][[time_idx]][temp_rthal_pt2[, 3] <= 0, ]
      }

      rthal_lpme_part_reconstructions <- reduce(
        rthal_lpme_part_reconstructions,
        rbind
      )

      rthal_pme_aug_reconstructions  <- reduce(
        rthal_pme_aug_reconstruction_list,
        rbind
      )

      rthal_pme_part_reconstructions <- map(
        seq_along(rthal_pme_part_reconstructions),
        ~ reduce(rthal_pme_part_reconstructions[[.x]], rbind)
      )
      rthal_pme_part_reconstructions <- reduce(
        rthal_pme_part_reconstructions, 
        rbind
      )

      rthal_pc_part_reconstructions <- map(
        seq_along(rthal_pc_part_reconstructions),
        ~ reduce(rthal_pc_part_reconstructions[[.x]], rbind)
      )
      rthal_pc_part_reconstructions <- reduce(
        rthal_pc_part_reconstructions, 
        rbind
      )

      rthal_lpme_aug_out <- list(
        lpme = rthal_lpme_aug,
        reconstructions = rthal_lpme_aug_reconstructions,
        fit_time = rthal_lpme_aug_time
      )
      rthal_pme_aug_out <- list(
        pme = rthal_pme_aug_list,
        reconstructions = rthal_pme_aug_reconstructions,
        fit_time = rthal_pme_aug_time_list
      )

      rthal_lpme_part_out <- list(
        lpme = rthal_lpme_part,
        reconstructions = rthal_lpme_part_reconstructions,
        fit_time = rthal_lpme_part_times
      )

      rthal_pme_part_out <- list(
        pme = rthal_pme_part,
        reconstructions = rthal_pme_part_reconstructions,
        fit_time = rthal_pme_part_times
      )

      rthal_pc_part_out <- list(
        pc = rthal_pc_part,
        reconstructions = rthal_pc_part_reconstructions,
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
    seed = TRUE
    # patno_idx = patno_idx,
    # patnos = patnos,
    # rthal_surface = rthal_surface
  )
}

rthal_results <- map(rthal_results, ~ value(.x))

