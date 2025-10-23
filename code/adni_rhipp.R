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

rhipp_surface <- read_csv(here("data/rhipp_surface_fsl.csv"))

adni_info <- read_csv(here("data/adni_info_full.csv")) |>
  distinct()
adni <- read_csv(here("data/adni_info.csv")) |>
  rename(
    image_id = "Image Data ID",
    subid = "Subject",
    date = "Acq Date"
  )

rhipp_processed <- preprocess_adni(rhipp_surface, adni_info)
rhipp_surface <- rhipp_processed$surface
rhipp_centers <- rhipp_processed$centers

write_csv(rhipp_surface, here("data/rhipp_surface_fsl_processed.csv"))

rhipp_fit <- fit_adni(
  rhipp_surface,
  rhipp_centers,
  "rhipp",
  verbose = FALSE,
  cores = parallel::detectCores() / 2
)
