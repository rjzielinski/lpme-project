fit_adni <- function(
  adni_surface,
  adni_centers,
  structure,
  partition_overlap = 0.15,
  cores = parallelly::availableCores(),
  verbose = FALSE
) {
  require(dplyr)
  require(future)
  require(future.batchtools)
  require(future.mirai)
  require(mirai)
  require(pme)
  require(tidyr)
  require(progressr)

  require(doRNG)
  require(doSNOW)
  require(foreach)

  source(here("code/functions/calculate_lpme_reconstructions.R"))
  source(here("code/functions/calculate_pme_reconstructions.R"))
  source(here("code/functions/estimate_volume.R"))
  source(here("code/functions/interior_identification.R"))
  source(here("code/prinSurf_v3.R"))
  source(here("code/functions/estimate_volume_interior.R"))

  handlers(global = TRUE)

  patnos <- unique(adni_surface$subid)
  patnos <- patnos[1:4]

  if (verbose == FALSE) {
    # plan(
    #   future.batchtools::batchtools_slurm,
    #   workers = 64,
    #   resources = list(
    #     time = "12:00:00",
    #     nodes = 1,
    #     ntasks = 1
    #   )
    # )
    # plan(cluster)

    # config <- cluster_config(
    #   command = "sbatch",
    #   options = "
    #   #SBATCH --job-name=adni_lpme_aug
    #   "
    # )
    #
    # daemons(n = cores, url = mirai::host_url(), remote = config)
    # while(mirai::info()[["connections"]] < cores) Sys.sleep(1.0)
    # daemons(cores)
    # plan(future.mirai::mirai_cluster)

    require(tictoc)
    require(pme)
    require(dplyr)
    require(tidyr)
    # daemons(cores)
    # everywhere(require(tictoc))
    # everywhere(require(pme))
    # everywhere(require(dplyr))
    # everywhere(require(tidyr))

    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  } else {
    plan(sequential, split = TRUE)
    require(tictoc)
    # daemons(1, output = TRUE)
    # everywhere(require(tictoc))
  }

  set.seed(500)

  adni_fit <- list()

  # p <- progressor(along = seq_along(patnos))
  pb <- txtProgressBar(max = length(patnos), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  registerDoRNG(716241)

  # for (patno_idx in seq_along(patnos)) {
  adni_fit <- foreach(
    patno_idx = seq_along(patnos),
    .export = c(
      "patnos",
      "adni_surface",
      "adni_centers",
      "structure",
      "partition_overlap",
      "verbose",
      "prinSurf",
      "calculate_lpme_reconstructions",
      "calculate_pme_reconstructions",
      "estimate_volume_interior_lpme",
      "estimate_volume_interior_pme"
    ),

    .packages = c(
      "alphashape3d",
      "dplyr",
      "here",
      "lubridate",
      "plotly",
      "plot3D",
      "pme",
      "pracma",
      "purrr",
      "RColorBrewer",
      "Rfast",
      "tictoc",
      "tidyr"
    ),
    .inorder = TRUE,
    .options.snow = opts,
    .errorhandling = "stop"
  ) %dopar%
    {
      patno <- patnos[patno_idx]
      patno_adni <- adni_surface |>
        filter(subid == patno) |>
        select(time_from_bl, x, y, z, theta, phi) |>
        mutate(
          partition = z > 0,
          partition1 = z > -partition_overlap,
          partition2 = z < partition_overlap
        )

      patno_adni_centers <- adni_centers |>
        filter(subid == patno) |>
        group_by(date) |>
        summarize(
          max_x = max(max_x),
          max_y = max(max_y),
          max_z = max(max_z)
        )

      adni_aug <- patno_adni |>
        select(-partition, -partition1, -partition2) |>
        as.matrix()

      adni_pt1 <- patno_adni |>
        filter(partition1 == TRUE) |>
        select(time_from_bl, x, y, z) |>
        as.matrix()
      adni_pt2 <- patno_adni |>
        filter(partition2 == TRUE) |>
        select(time_from_bl, x, y, z) |>
        as.matrix()

      time_values <- unique(patno_adni$time_from_bl)

      adni_pme_aug_list <- list()
      adni_pme_aug_time_list <- list()
      adni_pme_aug_reconstruction_list <- list()

      adni_lpme_part <- replicate(2, NULL, simplify = FALSE)
      adni_pme_part <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      adni_pc_part <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      adni_lpme_part_reconstructions <- replicate(
        2,
        NULL,
        simplify = FALSE
      )
      adni_pme_part_reconstructions <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      adni_pc_part_reconstructions <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      adni_lpme_part_times <- replicate(2, NULL, simplify = FALSE)
      adni_pme_part_times <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )
      adni_pc_part_times <- replicate(
        2,
        replicate(length(time_values), NULL, simplify = FALSE),
        simplify = FALSE
      )

      if (verbose) {
        print(paste0("Fitting augmented LPME"))
      }

      tic()
      adni_lpme_aug <- lpme(
        adni_aug,
        d = 2,
        verbose = FALSE,
        print_plots = FALSE
      )
      adni_lpme_aug_time <- toc(quiet = TRUE)
      adni_lpme_aug_time <- adni_lpme_aug_time$toc -
        adni_lpme_aug_time$tic

      adni_lpme_aug_reconstructions <- calculate_lpme_reconstructions(
        adni_lpme_aug,
        adni_aug
      )

      if (verbose) {
        print(paste0("Fitting partitioned LPME"))
      }

      tic()
      adni_lpme_part[[1]] <- lpme(
        adni_pt1,
        d = 2,
        verbose = FALSE,
        print_plots = FALSE
      )
      adni_lpme_time_pt1 <- toc(quiet = TRUE)
      adni_lpme_part_times[[1]] <- adni_lpme_time_pt1$toc -
        adni_lpme_time_pt1$tic
      adni_lpme_part_reconstructions[[
        1
      ]] <- calculate_lpme_reconstructions(
        adni_lpme_part[[1]],
        adni_pt1[adni_pt1[, 4] > 0, ]
      )
      tic()
      adni_lpme_part[[2]] <- lpme(
        adni_pt2,
        d = 2,
        verbose = FALSE,
        print_plots = FALSE
      )
      adni_lpme_time_pt2 <- toc(quiet = TRUE)
      adni_lpme_part_times[[2]] <- adni_lpme_time_pt2$toc -
        adni_lpme_time_pt2$tic
      adni_lpme_part_reconstructions[[
        2
      ]] <- calculate_lpme_reconstructions(
        adni_lpme_part[[2]],
        adni_pt2[adni_pt2[, 4] <= 0, ]
      )

      if (verbose) {
        print(paste0("Fitting time point-specific estimates"))
      }
      for (time_idx in seq_along(time_values)) {
        # print(time_idx)
        temp_adni <- adni_aug[adni_aug[, 1] == time_values[time_idx], -1]
        temp_adni_pt1 <- adni_pt1[
          adni_pt1[, 1] == time_values[time_idx],
          -1
        ]
        temp_adni_pt2 <- adni_pt2[
          adni_pt2[, 1] == time_values[time_idx],
          -1
        ]

        tic()
        adni_pme_aug_list[[time_idx]] <- pme(temp_adni, d = 2)
        temp_pme_aug_time <- toc(quiet = TRUE)
        adni_pme_aug_time_list[[time_idx]] <- temp_pme_aug_time$toc -
          temp_pme_aug_time$tic
        adni_pme_aug_reconstruction_list[[
          time_idx
        ]] <- calculate_pme_reconstructions(
          adni_pme_aug_list[[time_idx]],
          temp_adni
        )
        adni_pme_aug_reconstruction_list[[time_idx]] <- cbind(
          time_values[time_idx],
          adni_pme_aug_reconstruction_list[[time_idx]]
        )

        tic()
        adni_pme_part[[1]][[time_idx]] <- pme(temp_adni_pt1, d = 2)
        temp_pme_time_pt1 <- toc(quiet = TRUE)
        adni_pme_part_times[[1]][[time_idx]] <- temp_pme_time_pt1$toc -
          temp_pme_time_pt1$tic
        adni_pme_part_reconstructions[[1]][[
          time_idx
        ]] <- calculate_pme_reconstructions(
          adni_pme_part[[1]][[time_idx]],
          temp_adni_pt1[temp_adni_pt1[, 3] > 0, ]
        )
        tic()
        adni_pme_part[[2]][[time_idx]] <- pme(temp_adni_pt2, d = 2)
        temp_pme_time_pt2 <- toc(quiet = TRUE)
        adni_pme_part_times[[2]][[time_idx]] <- temp_pme_time_pt2$toc -
          temp_pme_time_pt2$tic
        adni_pme_part_reconstructions[[2]][[
          time_idx
        ]] <- calculate_pme_reconstructions(
          adni_pme_part[[2]][[time_idx]],
          temp_adni_pt2[temp_adni_pt2[, 3] <= 0, ]
        )

        adni_pme_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          adni_pme_part_reconstructions[[1]][[time_idx]]
        )
        adni_pme_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          adni_pme_part_reconstructions[[2]][[time_idx]]
        )

        tic()
        principal_surface_part1 <- prinSurf(temp_adni_pt1, flag = FALSE)
        temp_pc_time_pt1 <- toc(quiet = TRUE)
        adni_pc_part_times[[1]][[time_idx]] <- temp_pc_time_pt1$toc -
          temp_pc_time_pt1$tic
        surface_mse_part1 <- map(
          seq_along(principal_surface_part1),
          ~ principal_surface_part1[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part1 <- which.min(surface_mse_part1)

        adni_pc_part[[1]][[time_idx]] <- principal_surface_part1[[
          opt_surface_part1 + 2
        ]]
        adni_pc_part_reconstructions[[1]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part1[[opt_surface_part1 + 2]]$PS
        )
        adni_pc_part_reconstructions[[1]][[
          time_idx
        ]] <- adni_pc_part_reconstructions[[1]][[time_idx]][
          temp_adni_pt1[, 3] > 0,
        ]

        tic()
        principal_surface_part2 <- prinSurf(temp_adni_pt2, flag = FALSE)
        temp_pc_time_pt2 <- toc(quiet = TRUE)
        adni_pc_part_times[[2]][[time_idx]] <- temp_pc_time_pt2$toc -
          temp_pc_time_pt2$tic
        surface_mse_part2 <- map(
          seq_along(principal_surface_part2),
          ~ principal_surface_part2[[.x]]$MSE
        ) |>
          unlist()
        opt_surface_part2 <- which.min(surface_mse_part2)

        adni_pc_part[[2]][[time_idx]] <- principal_surface_part2[[
          opt_surface_part2 + 2
        ]]
        adni_pc_part_reconstructions[[2]][[time_idx]] <- cbind(
          time_values[time_idx],
          principal_surface_part2[[opt_surface_part2 + 2]]$PS
        )
        adni_pc_part_reconstructions[[2]][[
          time_idx
        ]] <- adni_pc_part_reconstructions[[2]][[time_idx]][
          temp_adni_pt2[, 3] <= 0,
        ]
      }

      if (verbose) {
        print(paste0("Fitting completed"))
      }

      adni_lpme_part_reconstructions <- reduce(
        adni_lpme_part_reconstructions,
        rbind
      )

      adni_lpme_part_times <- reduce(
        adni_lpme_part_times,
        sum
      )

      adni_pme_aug_reconstructions <- reduce(
        adni_pme_aug_reconstruction_list,
        rbind
      )

      adni_pme_aug_time_list <- reduce(
        adni_pme_aug_time_list,
        sum
      )

      adni_pme_part_reconstructions <- map(
        seq_along(adni_pme_part_reconstructions),
        ~ reduce(adni_pme_part_reconstructions[[.x]], rbind)
      )
      adni_pme_part_reconstructions <- reduce(
        adni_pme_part_reconstructions,
        rbind
      )

      adni_pme_part_times <- map(
        seq_along(adni_pme_part_times),
        ~ reduce(adni_pme_part_times[[.x]], sum)
      )
      adni_pme_part_times <- reduce(
        adni_pme_part_times,
        sum
      )

      adni_pc_part_reconstructions <- map(
        seq_along(adni_pc_part_reconstructions),
        ~ reduce(adni_pc_part_reconstructions[[.x]], rbind)
      )
      adni_pc_part_reconstructions <- reduce(
        adni_pc_part_reconstructions,
        rbind
      )

      adni_pc_part_times <- map(
        seq_along(adni_pc_part_times),
        ~ reduce(adni_pc_part_times[[.x]], sum)
      )
      adni_pc_part_times <- reduce(
        adni_pc_part_times,
        sum
      )

      if (verbose) {
        print(paste0("Estimating volumes"))
      }

      adni_lpme_part_volumes <- estimate_volume_interior_lpme(
        adni_lpme_part,
        list(adni_pt1, adni_pt2),
        time_values,
        n_points = 10000,
        data_max = patno_adni_centers,
        limit_scaler = 0.05,
        partition_index = 3
      )
      adni_lpme_part_volumes <- adni_lpme_part_volumes$volumes

      adni_pme_part_volumes <- estimate_volume_interior_pme(
        adni_pme_part,
        list(adni_pt1, adni_pt2),
        time_values,
        n_points = 10000,
        data_max = patno_adni_centers,
        limit_scaler = 0.05,
        partition_index = 3
      )
      adni_pme_part_volumes <- adni_pme_part_volumes$volumes

      adni_lpme_aug_volumes <- vector()
      adni_pme_aug_volumes <- vector()
      adni_pc_part_volumes <- vector()

      adni_lpme_aug_interior <- list()
      adni_pme_aug_interior <- list()
      adni_pc_part_interior <- list()

      adni_lpme_part_volumes_mesh <- vector()
      adni_pme_part_volumes_mesh <- vector()

      for (time_idx in seq_along(time_values)) {
        x_scale <- patno_adni_centers$max_x[time_idx]
        y_scale <- patno_adni_centers$max_y[time_idx]
        z_scale <- patno_adni_centers$max_z[time_idx]
        temp_lpme_reconstructions <- adni_lpme_aug_reconstructions[
          adni_lpme_aug_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pme_reconstructions <- adni_pme_aug_reconstructions[
          adni_pme_aug_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pc_reconstructions <- adni_pc_part_reconstructions[
          adni_pc_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]

        temp_lpme_part_reconstructions <- adni_lpme_part_reconstructions[
          adni_lpme_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]
        temp_pme_part_reconstructions <- adni_pme_part_reconstructions[
          adni_pme_part_reconstructions[, 1] == time_values[time_idx],
          2:4
        ]

        temp_lpme_ashape <- ashape3d(
          temp_lpme_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        adni_lpme_aug_volumes[time_idx] <- volume_ashape3d(
          temp_lpme_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pme_ashape <- ashape3d(
          temp_pme_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        adni_pme_aug_volumes[time_idx] <- volume_ashape3d(
          temp_pme_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pc_ashape <- ashape3d(
          temp_pc_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        adni_pc_part_volumes[time_idx] <- volume_ashape3d(
          temp_pc_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_lpme_part_ashape <- ashape3d(
          temp_lpme_part_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        adni_lpme_part_volumes_mesh[time_idx] <- volume_ashape3d(
          temp_lpme_part_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)

        temp_pme_part_ashape <- ashape3d(
          temp_pme_part_reconstructions,
          alpha = seq(0.5, 1.5, by = 0.1)
        )
        adni_pme_part_volumes_mesh[time_idx] <- volume_ashape3d(
          temp_pme_part_ashape,
          indexAlpha = "all"
        ) |>
          median() *
          (x_scale * y_scale * z_scale)
      }

      adni_lpme_aug_out <- list(
        lpme = adni_lpme_aug,
        reconstructions = adni_lpme_aug_reconstructions,
        volumes = adni_lpme_aug_volumes,
        fit_time = adni_lpme_aug_time
      )
      adni_pme_aug_out <- list(
        pme = adni_pme_aug_list,
        reconstructions = adni_pme_aug_reconstructions,
        volumes = adni_pme_aug_volumes,
        fit_time = adni_pme_aug_time_list
      )

      adni_lpme_part_out <- list(
        lpme = adni_lpme_part,
        reconstructions = adni_lpme_part_reconstructions,
        volumes = adni_lpme_part_volumes,
        volumes_mesh = adni_lpme_part_volumes_mesh,
        fit_time = adni_lpme_part_times
      )

      adni_pme_part_out <- list(
        pme = adni_pme_part,
        reconstructions = adni_pme_part_reconstructions,
        volumes = adni_pme_part_volumes,
        volumes_mesh = adni_pme_part_volumes_mesh,
        fit_time = adni_pme_part_times
      )

      adni_pc_part_out <- list(
        pc = adni_pc_part,
        reconstructions = adni_pc_part_reconstructions,
        volumes = adni_pc_part_volumes,
        fit_time = adni_pc_part_times
      )

      adni <- list(
        lpme_aug = adni_lpme_aug_out,
        pme_aug = adni_pme_aug_out,
        lpme_part = adni_lpme_part_out,
        pme_part = adni_pme_part_out,
        pc_part = adni_pc_part_out
      )

      if (!dir.exists(here(paste0("output/adni/", patno)))) {
        dir.create(here(paste0("output/adni/", patno)), recursive = TRUE)
      }

      print(
        paste0(
          "Saving results for participant ",
          patno,
          ", structure ",
          structure
        )
      )

      saveRDS(
        adni,
        here(paste0(
          "output/adni/",
          patno,
          "/",
          structure,
          "_results.RDS"
        ))
      )
      TRUE
      # p(sprintf("subject number: %g", patno_idx))
    }
  # seed = TRUE
  # patno_adni = patno_adni,
  # patno_adni_centers = patno_adni_centers,
  # adni_aug = adni_aug,
  # adni_pt1 = adni_pt1,
  # adni_pt2 = adni_pt2,
  # time_values = time_values,
  # verbose = verbose
  #   )
  # }
  #
  # # collect_mirai(adni_fit)
  # results <- map(adni_fit, value)

  results <- adni_fit

  # plan(sequential)
  # invisible(daemons(0))
}
