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

lhipp_surface <- read_csv(here("data/lhipp_surface_fsl.csv"))

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

write_csv(lhipp_surface, here("data/lhipp_surface_fsl_processed.csv"))


lhipp_fit <- fit_adni(
  lhipp_surface,
  lhipp_centers,
  "lhipp",
  verbose = FALSE,
  cores = parallel::detectCores() / 2
)
