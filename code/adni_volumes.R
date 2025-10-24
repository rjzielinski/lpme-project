library(dplyr)
library(foreach)
library(here)
library(plotly)
library(pme)
library(readr)
library(reticulate)
library(RColorBrewer)

library(doRNG)
library(doSNOW)

use_condaenv("lpme")
np <- import("numpy")
o3d <- import("open3d")
pv <- import("pyvista")

# source(here("code/functions/estimate_mesh_volume.R"))
source(here("code/functions/estimate_mesh_volume_poisson.R"))

cores <- parallel::detectCores() - 2

result_dir <- here("output/adni")

isolate_patno <- function(file_vec, file_path, structure) {
  patnos <- gsub(file_path, "", file_vec)
  patnos <- gsub(paste0(structure, "_results.RDS"), "", patnos)
  patnos <- gsub("/", "", patnos)
}

result_files <- list.files(result_dir, full.names = TRUE, recursive = TRUE)

lhipp_files <- result_files[grepl("lhipp", result_files)]
rhipp_files <- result_files[grepl("rhipp", result_files)]
lthal_files <- result_files[grepl("lthal", result_files)]
rthal_files <- result_files[grepl("rthal", result_files)]

lhipp_patnos <- isolate_patno(lhipp_files, result_dir, "lhipp")
rhipp_patnos <- isolate_patno(rhipp_files, result_dir, "rhipp")
lthal_patnos <- isolate_patno(lthal_files, result_dir, "lthal")
rthal_patnos <- isolate_patno(rthal_files, result_dir, "rthal")

common_patnos <- intersect(lhipp_patnos, lthal_patnos)
common_patnos <- intersect(common_patnos, rhipp_patnos)
common_patnos <- intersect(common_patnos, rthal_patnos)

# RETRIEVE VOLUME ESTIMATES FROM MESHES -------------------------------------

adni <- read_csv("data/adni_info_full.csv")

lhipp_surface <- read_csv(here("data/lhipp_surface_fsl_processed.csv"))
rhipp_surface <- read_csv(here("data/rhipp_surface_fsl_processed.csv"))
lthal_surface <- read_csv(here("data/lthal_surface_fsl_processed.csv"))
rthal_surface <- read_csv(here("data/rthal_surface_fsl_processed.csv"))

cl <- makeCluster(cores, outfile = "log.txt")
registerDoSNOW(cl)
registerDoRNG(42)

adni_volumes <- foreach(
  patno = common_patnos,
  .combine = "rbind",
  .inorder = FALSE,
  .export = c(
    "common_patnos",
    "adni",
    "lhipp_surface",
    "rhipp_surface",
    "lthal_surface",
    "rthal_surface",
    "result_dir",
    "estimate_mesh_volume_poisson"
  ),
  .packages = c(
    "dplyr",
    "here",
    "pme",
    "purrr",
    "reticulate",
    "tidyr"
  )
) %dopar%
  {
    use_condaenv("lpme")
    np <- import("numpy")
    o3d <- import("open3d")
    pv <- import("pyvista")

    print(
      paste0(
        "PATNO ",
        which(common_patnos == patno),
        " of ",
        length(common_patnos),
        ": ",
        patno
      )
    )

    patno_dir <- paste0(result_dir, "/", patno, "/")
    lhipp_info <- lhipp_surface |>
      filter(subid == patno)
    rhipp_info <- rhipp_surface |>
      filter(subid == patno)
    lthal_info <- lthal_surface |>
      filter(subid == patno)
    rthal_info <- rthal_surface |>
      filter(subid == patno)

    gc()

    patno_dates <- unique(c(
      lhipp_info$date,
      rhipp_info$date,
      lthal_info$date,
      rthal_info$date
    ))

    lhipp_result <- readRDS(paste0(patno_dir, "lhipp_results.RDS"))
    rhipp_result <- readRDS(paste0(patno_dir, "rhipp_results.RDS"))
    lthal_result <- readRDS(paste0(patno_dir, "lthal_results.RDS"))
    rthal_result <- readRDS(paste0(patno_dir, "rthal_results.RDS"))

    time_points <- unique(lhipp_result$data$time_from_bl)
    lhipp_data_vol <- vector(mode = "numeric", length = length(time_points))
    rhipp_data_vol <- vector(mode = "numeric", length = length(time_points))
    lthal_data_vol <- vector(mode = "numeric", length = length(time_points))
    rthal_data_vol <- vector(mode = "numeric", length = length(time_points))

    for (time_idx in seq_along(time_points)) {
      lhipp_scaling <- lhipp_info |>
        filter(time_from_bl == time_points[time_idx]) |>
        summarize(max_x = mean(max_x), max_y = mean(max_y), max_z = mean(max_z))
      rhipp_scaling <- rhipp_info |>
        filter(time_from_bl == time_points[time_idx]) |>
        summarize(max_x = mean(max_x), max_y = mean(max_y), max_z = mean(max_z))
      lthal_scaling <- lthal_info |>
        filter(time_from_bl == time_points[time_idx]) |>
        summarize(max_x = mean(max_x), max_y = mean(max_y), max_z = mean(max_z))
      rthal_scaling <- rthal_info |>
        filter(time_from_bl == time_points[time_idx]) |>
        summarize(max_x = mean(max_x), max_y = mean(max_y), max_z = mean(max_z))

      temp_lhipp_data <- lhipp_result$data |>
        filter(time_from_bl == time_points[time_idx]) |>
        select(x, y, z) |>
        as.matrix()
      temp_rhipp_data <- rhipp_result$data |>
        filter(time_from_bl == time_points[time_idx]) |>
        select(x, y, z) |>
        as.matrix()

      temp_lthal_data <- lthal_result$data |>
        filter(time_from_bl == time_points[time_idx]) |>
        select(x, y, z) |>
        as.matrix()
      temp_rthal_data <- rthal_result$data |>
        filter(time_from_bl == time_points[time_idx]) |>
        select(x, y, z) |>
        as.matrix()

      lhipp_data_vol[time_idx] <- estimate_mesh_volume_poisson(
        temp_lhipp_data
      )$volume *
        (lhipp_scaling$max_x *
          lhipp_scaling$max_y *
          lhipp_scaling$max_z)
      rhipp_data_vol[time_idx] <- estimate_mesh_volume_poisson(
        temp_rhipp_data
      )$volume *
        (rhipp_scaling$max_x *
          rhipp_scaling$max_y *
          rhipp_scaling$max_z)

      lthal_data_vol[time_idx] <- estimate_mesh_volume_poisson(
        temp_lthal_data,
        alpha_vals = exp(seq(-2, 0, 0.25))
      )$volume *
        (lthal_scaling$max_x *
          lthal_scaling$max_y *
          lthal_scaling$max_z)
      rthal_data_vol[time_idx] <- estimate_mesh_volume_poisson(
        temp_rthal_data,
        alpha_vals = exp(seq(-2, 0, 0.25))
      )$volume *
        (rthal_scaling$max_x *
          rthal_scaling$max_y *
          rthal_scaling$max_z)
    }

    lhipp_lpme_vols <- lhipp_result$lpme_part$volumes_mesh
    lhipp_pme_vols <- lhipp_result$pme_part$volumes_mesh
    lhipp_pc_vols <- lhipp_result$pc_part$volumes

    rhipp_lpme_vols <- rhipp_result$lpme_part$volumes_mesh
    rhipp_pme_vols <- rhipp_result$pme_part$volumes_mesh
    rhipp_pc_vols <- rhipp_result$pc_part$volumes

    lthal_lpme_vols <- lthal_result$lpme_part$volumes_mesh
    lthal_pme_vols <- lthal_result$pme_part$volumes_mesh
    lthal_pc_vols <- lthal_result$pc_part$volumes

    rthal_lpme_vols <- rthal_result$lpme_part$volumes_mesh
    rthal_pme_vols <- rthal_result$pme_part$volumes_mesh
    rthal_pc_vols <- rthal_result$pc_part$volumes

    patno_info <- data.frame(
      patno = rep(patno, length(patno_dates)),
      date = patno_dates,
      lhipp_data_vol = lhipp_data_vol,
      lhipp_lpme_vol = lhipp_lpme_vols,
      lhipp_pme_vol = lhipp_pme_vols,
      lhipp_pc_vol = lhipp_pc_vols,
      rhipp_data_vol = rhipp_data_vol,
      rhipp_lpme_vol = rhipp_lpme_vols,
      rhipp_pme_vol = rhipp_pme_vols,
      rhipp_pc_vol = rhipp_pc_vols,
      lthal_data_vol = lthal_data_vol,
      lthal_lpme_vol = lthal_lpme_vols,
      lthal_pme_vol = lthal_pme_vols,
      lthal_pc_vol = lthal_pc_vols,
      rthal_data_vol = rthal_data_vol,
      rthal_lpme_vol = rthal_lpme_vols,
      rthal_pme_vol = rthal_pme_vols,
      rthal_pc_vol = rthal_pc_vols
    )
  }

stopCluster(cl)

write_csv(adni_volumes, here("output/adni_volumes.csv"))
