# library(alphashape3d)
library(here)
library(lubridate)
# library(mirai)
library(future)
library(plotly)
# library(plot3D)
library(pme)
library(pracma)
library(RColorBrewer)
library(Rfast)
library(tictoc)
library(tidyverse)

library(reticulate)
use_condaenv("lpme")

source(here("code/functions/preprocess_adni.R"))
source(here("code/functions/fit_adni.R"))

# options(future.globals.maxSize = 4 * 1e9)

voxel_vol <- 1.2 * 0.9375 * 0.9375

rthal_surface <- read_csv(here("data/rthal_surface_fsl.csv"))

adni_info <- read_csv(here("data/adni_info_full.csv")) |>
  distinct()
adni <- read_csv(here("data/adni_info.csv")) |>
  rename(
    image_id = "Image Data ID",
    subid = "Subject",
    date = "Acq Date"
  )

rthal_processed <- preprocess_adni(rthal_surface, adni_info)
rthal_surface <- rthal_processed$surface
rthal_centers <- rthal_processed$centers

write_csv(rthal_surface, here("data/rthal_surface_fsl_processed.csv"))

rthal_fit <- fit_adni(
  rthal_surface,
  rthal_centers,
  "rthal",
  verbose = FALSE,
  cores = parallel::detectCores() / 2
  # verbose = TRUE
)
