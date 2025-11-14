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

# cores <- parallel::detectCores() - 2

result_dir <- here("output/adni_area")

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

common_patnos <- intersect(lhipp_patnos, lthal_patnos) |>
  intersect(rhipp_patnos) |>
  intersect(rthal_patnos)

# RETRIEVE VOLUME ESTIMATES FROM MESHES -------------------------------------

adni <- read_csv("data/adni_info_full.csv")

lhipp_surface <- read_csv(
  here("data/lhipp_surface_fsl_processed.csv"),
  col_select = c("subid", "date")
) |>
  distinct()
gc()

rhipp_surface <- read_csv(
  here("data/rhipp_surface_fsl_processed.csv"),
  col_select = c("subid", "date")
) |>
  distinct()
gc()

lthal_surface <- read_csv(
  here("data/lthal_surface_fsl_processed.csv"),
  col_select = c("subid", "date")
) |>
  distinct()
gc()

rthal_surface <- read_csv(
  here("data/rthal_surface_fsl_processed.csv"),
  col_select = c("subid", date)
) |>
  distinct()
gc()

# cl <- makeCluster(cores, outfile = "log.txt")
# registerDoSNOW(cl)
# registerDoRNG(42)

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
  ),
  .errorhandling = "remove"
) %do%
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

    time_points <- unique(lhipp_result$data$data$time_from_bl)
    lhipp_data_area_dim1 <- lhipp_result$data$area_dim1
    lhipp_data_area_dim2 <- lhipp_result$data$area_dim2
    lhipp_data_area_dim3 <- lhipp_result$data$area_dim3

    rhipp_data_area_dim1 <- rhipp_result$data$area_dim1
    rhipp_data_area_dim2 <- rhipp_result$data$area_dim2
    rhipp_data_area_dim3 <- rhipp_result$data$area_dim3

    lthal_data_area_dim1 <- lthal_result$data$area_dim1
    lthal_data_area_dim2 <- lthal_result$data$area_dim2
    lthal_data_area_dim3 <- lthal_result$data$area_dim3

    rthal_data_area_dim1 <- rthal_result$data$area_dim1
    rthal_data_area_dim2 <- rthal_result$data$area_dim2
    rthal_data_area_dim3 <- rthal_result$data$area_dim3

    lhipp_lpme_part_area_dim1 <- lhipp_result$lpme_part$area_dim1
    lhipp_lpme_part_area_dim2 <- lhipp_result$lpme_part$area_dim2
    lhipp_lpme_part_area_dim3 <- lhipp_result$lpme_part$area_dim3
    lhipp_lpme_part_time <- lhipp_result$lpme_part$fit_time

    lhipp_pme_part_area_dim1 <- lhipp_result$pme_part$area_dim1
    lhipp_pme_part_area_dim2 <- lhipp_result$pme_part$area_dim2
    lhipp_pme_part_area_dim3 <- lhipp_result$pme_part$area_dim3
    lhipp_pme_part_time <- lhipp_result$pme_part$fit_time

    lhipp_pc_part_area_dim1 <- lhipp_result$pc_part$area_dim1
    lhipp_pc_part_area_dim2 <- lhipp_result$pc_part$area_dim2
    lhipp_pc_part_area_dim3 <- lhipp_result$pc_part$area_dim3
    lhipp_pc_part_time <- lhipp_result$pc_part$fit_time

    rhipp_lpme_part_area_dim1 <- rhipp_result$lpme_part$area_dim1
    rhipp_lpme_part_area_dim2 <- rhipp_result$lpme_part$area_dim2
    rhipp_lpme_part_area_dim3 <- rhipp_result$lpme_part$area_dim3
    rhipp_lpme_part_time <- rhipp_result$lpme_part$fit_time

    rhipp_pme_part_area_dim1 <- rhipp_result$pme_part$area_dim1
    rhipp_pme_part_area_dim2 <- rhipp_result$pme_part$area_dim2
    rhipp_pme_part_area_dim3 <- rhipp_result$pme_part$area_dim3
    rhipp_pme_part_time <- rhipp_result$pme_part$fit_time

    rhipp_pc_part_area_dim1 <- rhipp_result$pc_part$area_dim1
    rhipp_pc_part_area_dim2 <- rhipp_result$pc_part$area_dim2
    rhipp_pc_part_area_dim3 <- rhipp_result$pc_part$area_dim3
    rhipp_pc_part_time <- rhipp_result$pc_part$fit_time

    lthal_lpme_part_area_dim1 <- lthal_result$lpme_part$area_dim1
    lthal_lpme_part_area_dim2 <- lthal_result$lpme_part$area_dim2
    lthal_lpme_part_area_dim3 <- lthal_result$lpme_part$area_dim3
    lthal_lpme_part_time <- lthal_result$lpme_part$fit_time

    lthal_pme_part_area_dim1 <- lthal_result$pme_part$area_dim1
    lthal_pme_part_area_dim2 <- lthal_result$pme_part$area_dim2
    lthal_pme_part_area_dim3 <- lthal_result$pme_part$area_dim3
    lthal_pme_part_time <- lthal_result$pme_part$fit_time

    lthal_pc_part_area_dim1 <- lthal_result$pc_part$area_dim1
    lthal_pc_part_area_dim2 <- lthal_result$pc_part$area_dim2
    lthal_pc_part_area_dim3 <- lthal_result$pc_part$area_dim3
    lthal_pc_part_time <- lthal_result$pc_part$fit_time

    rthal_lpme_part_area_dim1 <- rthal_result$lpme_part$area_dim1
    rthal_lpme_part_area_dim2 <- rthal_result$lpme_part$area_dim2
    rthal_lpme_part_area_dim3 <- rthal_result$lpme_part$area_dim3
    rthal_lpme_part_time <- rthal_result$lpme_part$fit_time

    rthal_pme_part_area_dim1 <- rthal_result$pme_part$area_dim1
    rthal_pme_part_area_dim2 <- rthal_result$pme_part$area_dim2
    rthal_pme_part_area_dim3 <- rthal_result$pme_part$area_dim3
    rthal_pme_part_time <- rthal_result$pme_part$fit_time

    rthal_pc_part_area_dim1 <- rthal_result$pc_part$area_dim1
    rthal_pc_part_area_dim2 <- rthal_result$pc_part$area_dim2
    rthal_pc_part_area_dim3 <- rthal_result$pc_part$area_dim3
    rthal_pc_part_time <- rthal_result$pc_part$fit_time

    patno_info <- data.frame(
      patno = rep(patno, length(patno_dates)),
      date = patno_dates,
      lhipp_data_area_dim1 = lhipp_data_area_dim1,
      lhipp_data_area_dim2 = lhipp_data_area_dim2,
      lhipp_data_area_dim3 = lhipp_data_area_dim3,
      lhipp_lpme_area_dim1 = lhipp_lpme_part_area_dim1,
      lhipp_lpme_area_dim2 = lhipp_lpme_part_area_dim2,
      lhipp_lpme_area_dim3 = lhipp_lpme_part_area_dim3,
      lhipp_lpme_time = lhipp_lpme_part_time,
      lhipp_pme_area_dim1 = lhipp_pme_part_area_dim1,
      lhipp_pme_area_dim2 = lhipp_pme_part_area_dim2,
      lhipp_pme_area_dim3 = lhipp_pme_part_area_dim3,
      lhipp_pme_time = lhipp_pme_part_time,
      lhipp_pc_area_dim1 = lhipp_pc_part_area_dim1,
      lhipp_pc_area_dim2 = lhipp_pc_part_area_dim2,
      lhipp_pc_area_dim3 = lhipp_pc_part_area_dim3,
      lhipp_pc_time = lhipp_pc_part_time,
      rhipp_data_area_dim1 = rhipp_data_area_dim1,
      rhipp_data_area_dim2 = rhipp_data_area_dim2,
      rhipp_data_area_dim3 = rhipp_data_area_dim3,
      rhipp_lpme_area_dim1 = rhipp_lpme_part_area_dim1,
      rhipp_lpme_area_dim2 = rhipp_lpme_part_area_dim2,
      rhipp_lpme_area_dim3 = rhipp_lpme_part_area_dim3,
      rhipp_lpme_time = rhipp_lpme_part_time,
      rhipp_pme_area_dim1 = rhipp_pme_part_area_dim1,
      rhipp_pme_area_dim2 = rhipp_pme_part_area_dim2,
      rhipp_pme_area_dim3 = rhipp_pme_part_area_dim3,
      rhipp_pme_time = rhipp_pme_part_time,
      rhipp_pc_area_dim1 = rhipp_pc_part_area_dim1,
      rhipp_pc_area_dim2 = rhipp_pc_part_area_dim2,
      rhipp_pc_area_dim3 = rhipp_pc_part_area_dim3,
      rhipp_pc_time = rhipp_pc_part_time,
      lthal_data_area_dim1 = lthal_data_area_dim1,
      lthal_data_area_dim2 = lthal_data_area_dim2,
      lthal_data_area_dim3 = lthal_data_area_dim3,
      lthal_lpme_area_dim1 = lthal_lpme_part_area_dim1,
      lthal_lpme_area_dim2 = lthal_lpme_part_area_dim2,
      lthal_lpme_area_dim3 = lthal_lpme_part_area_dim3,
      lthal_lpme_time = lthal_lpme_part_time,
      lthal_pme_area_dim1 = lthal_pme_part_area_dim1,
      lthal_pme_area_dim2 = lthal_pme_part_area_dim2,
      lthal_pme_area_dim3 = lthal_pme_part_area_dim3,
      lthal_pme_time = lthal_pme_part_time,
      lthal_pc_area_dim1 = lthal_pc_part_area_dim1,
      lthal_pc_area_dim2 = lthal_pc_part_area_dim2,
      lthal_pc_area_dim3 = lthal_pc_part_area_dim3,
      lthal_pc_time = lthal_pc_part_time,
      rthal_data_area_dim1 = rthal_data_area_dim1,
      rthal_data_area_dim2 = rthal_data_area_dim2,
      rthal_data_area_dim3 = rthal_data_area_dim3,
      rthal_lpme_area_dim1 = rthal_lpme_part_area_dim1,
      rthal_lpme_area_dim2 = rthal_lpme_part_area_dim2,
      rthal_lpme_area_dim3 = rthal_lpme_part_area_dim3,
      rthal_lpme_time = rthal_lpme_part_time,
      rthal_pme_area_dim1 = rthal_pme_part_area_dim1,
      rthal_pme_area_dim2 = rthal_pme_part_area_dim2,
      rthal_pme_area_dim3 = rthal_pme_part_area_dim3,
      rthal_pme_time = rthal_pme_part_time,
      rthal_pc_area_dim1 = rthal_pc_part_area_dim1,
      rthal_pc_area_dim2 = rthal_pc_part_area_dim2,
      rthal_pc_area_dim3 = rthal_pc_part_area_dim3,
      rthal_pc_time = rthal_pc_part_time
    )
  }

# stopCluster(cl)

write_csv(adni_volumes, here("output/adni_areas.csv"))
