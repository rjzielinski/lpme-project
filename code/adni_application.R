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
  lhipp_surface,
  lhipp_centers,
  "lhipp",
  cores = 6,
  verbose = FALSE
)
daemons(0)
rhipp_fit <- fit_adni(
  rhipp_surface,
  rhipp_centers,
  "rhipp",
  cores = 6,
  verbose = FALSE
)
daemons(0)
lthal_fit <- fit_adni(
  lthal_surface,
  lthal_centers,
  "lthal",
  cores = 6,
  verbose = FALSE
)
daemons(0)
rthal_fit <- fit_adni(
  rthal_surface,
  rthal_centers,
  "rthal",
  cores = 6,
  verbose = FALSE
)
daemons(0)
