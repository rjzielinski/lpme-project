library(oro.dicom)
library(oro.nifti)
library(neurobase)
library(fslr)
library(extrantsr)
library(stringr)
library(foreach)
library(tidyverse)

sub_dirs <- list.dirs("data/ADNI", recursive = FALSE)

nii_scans <- list()

processed_scans <- list()

hipp_mask <- readNIfTI(
  paste0(fsl_dir(), "/data/standard/MNI152_T1_1mm_Hipp_mask_dil8.nii.gz")
)

ncores <- parallel::detectCores() / 2
cl <- parallel::makeCluster(ncores, type = "FORK")
doParallel::registerDoParallel(cl = cl)

foreach (dir_idx = 1:length(sub_dirs)) %do% {
  img_dirs <- list.files(sub_dirs[dir_idx], recursive = TRUE, full.names = TRUE) %>%
    gsub(pattern = "([^/]+$)", replacement = "") %>%
    unique()

  nii_scans[[dir_idx]] <- list()

  for (img_idx in 1:length(img_dirs)) {
    all_slices <- readDICOM(img_dirs[img_idx])
    nii <- dicom2nifti(all_slices)

    # Inhomogeneity correction
    # Brain extraction / segmentation
    # Image Registration
    # Tissue Class Segmentation

    nii_scans[[dir_idx]][[img_idx]] <- nii

    # inhomogeneity correction
    # nii_bc <- bias_correct(nii, correction = "N4")

    # brain extraction
    # nii_bc_be <- fslbet(nii_bc)
    # if brain extraction shows signs of issues in some images,
    # extrantsr::fslbet_robust() may be useful
    # includes options for neck correction, etc.

    # register image to MNI space
    # not using coregistration to allow for future population-wide analysis
    # nii_reg <- registration(
    #   nii,
    #   skull_strip = FALSE,
    #   correct = FALSE
    # )

    proc_dir <- paste0("data/ADNI_processed", gsub(pattern = "data/ADNI", replacement = "", img_dirs[img_idx]))

    run_first_all(
      # nii_reg$outfile,
      nii,
      oprefix = proc_dir,
      brain_extracted = TRUE
    )
  }
}

parallel::stopCluster(cl = cl)

test_dicom <- readDICOM(
  list.files(sub_dirs[1], recursive = TRUE, full.names = TRUE)[1]
)

patnos <- list.dirs("data/ADNI_processed", recursive = FALSE, full.names = FALSE)
processed_dirs <- list.dirs("data/ADNI_processed", recursive = FALSE)

lhipp <- matrix(ncol = 6)
rhipp <- matrix(ncol = 6)

for (dir_idx in 1:length(processed_dirs)) {
  print(patnos[dir_idx])
  scan_dates <- list.dirs(
    paste0(processed_dirs[dir_idx], "/MP-RAGE"),
    full.names = FALSE,
    recursive = FALSE
  )
  scan_dirs <- list.files(
    paste0(processed_dirs[dir_idx], "/MP-RAGE"),
    recursive = TRUE,
    full.names = TRUE
  ) %>%
    gsub(pattern = "([^/]+$)", replacement = "") %>%
    unique()
  for (scan_idx in 1:length(scan_dirs)) {
    print(scan_dates[scan_idx])
    temp_lhipp <- readnii(paste0(scan_dirs[scan_idx], "-L_Hipp_first.nii.gz"))
    temp_rhipp <- readnii(paste0(scan_dirs[scan_idx], "-R_Hipp_first.nii.gz"))

    lhipp_idx <- which(temp_lhipp > 0, arr.ind = TRUE)
    rhipp_idx <- which(temp_rhipp > 0, arr.ind = TRUE)

    lhipp_idx <- cbind(patnos[dir_idx], scan_dates[scan_idx], lhipp_idx, temp_lhipp[lhipp_idx])
    rhipp_idx <- cbind(patnos[dir_idx], scan_dates[scan_idx], rhipp_idx, temp_rhipp[rhipp_idx])

    lhipp <- rbind(lhipp, lhipp_idx)
    rhipp <- rbind(rhipp, rhipp_idx)
  }
}

lhipp <- lhipp[-1, ] %>%
  as_tibble() %>%
  rename(
    patno = V1,
    scan_date = V2,
    x = dim1,
    y = dim2,
    z = dim3,
    intensity = V6
  ) %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y),
    z = as.numeric(z),
    intensity = as.numeric(intensity)
  )
rhipp <- rhipp[-1, ] %>%
  as_tibble() %>%
  rename(
    patno = V1,
    scan_date = V2,
    x = dim1,
    y = dim2,
    z = dim3,
    intensity = V6
  ) %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y),
    z = as.numeric(z),
    intensity = as.numeric(intensity)
  )

write_csv(lhipp, "data/adni_lhipp.csv")
write_csv(rhipp, "data/adni_rhipp.csv")