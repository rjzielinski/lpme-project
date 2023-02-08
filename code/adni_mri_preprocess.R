library(oro.dicom)
library(oro.nifti)
library(neurobase)
library(fslr)
library(extrantsr)
library(stringr)
library(foreach)
library(tidyverse)

sub_dirs <- list.dirs("data/adni", recursive = FALSE)

nii_scans <- list()

processed_scans <- list()

hipp_mask <- readNIfTI(
  paste0(fsl_dir(), "/data/standard/MNI152_T1_1mm_Hipp_mask_dil8.nii.gz")
)

mni_nifti <- readNIfTI(
  paste0(fsl_dir(), "/data/standard/MNI152_T1_1mm.nii.gz")
)

ncores <- parallel::detectCores() / 2
cl <- parallel::makeCluster(ncores, type = "FORK")
doParallel::registerDoParallel(cl = cl)

foreach (dir_idx = 1:length(sub_dirs)) %dopar% {
  img_dirs <- list.files(sub_dirs[dir_idx], recursive = TRUE, full.names = TRUE) %>%
    gsub(pattern = "([^/]+$)", replacement = "") %>%
    unique()

  nii_scans[[dir_idx]] <- list()

  for (img_idx in 1:length(img_dirs)) {
    proc_dir <- paste0("data/adni_processed_fsl", gsub(pattern = "data/adni", replacement = "", img_dirs[img_idx]))
    if (file.exists(paste0(proc_dir, "/-L_Hipp_first.nii.gz"))) {
      next
    } else {
      all_slices <- readDICOM(img_dirs[img_idx])
      nii <- dicom2nifti(all_slices)

      # Inhomogeneity correction
      # Brain extraction / segmentation
      # Image Registration
      # Tissue Class Segmentation

      # nii_scans[[dir_idx]][[img_idx]] <- nii

      nii <- fslreorient2std(nii, verbose = FALSE) # orient nii to mni space
      # inhomogeneity correction
      nii <- bias_correct(
        nii,
        correction = "N4",
        verbose = FALSE
      ) # N4 bias correction using ANTs
      # nii_bc <- fast(nii, bias_correct = TRUE)

      # brain extraction
      # if brain extraction shows signs of issues in some images,
      # extrantsr::fslbet_robust() may be useful
      # includes options for neck correction, etc.

      # register image to MNI space
      # not using coregistration to allow for future population-wide analysis

      # nii_reg <- registration(
      #   nii,
      #   skull_strip = FALSE,
      #   correct = TRUE
      # )

      # Do not run brain extraction -- recommended that first_flirt runs on full-head images
      # nii <- fslbet(nii, retimg = TRUE)

      run_first_all(
        nii,
        oprefix = proc_dir,
        verbose = FALSE
      )
    }
  }
}

parallel::stopCluster(cl = cl)

test_dicom <- readDICOM(
  list.files(sub_dirs[1], recursive = TRUE, full.names = TRUE)[1]
)

patnos <- list.dirs("data/adni_processed_fsl", recursive = FALSE, full.names = FALSE)
processed_dirs <- list.dirs("data/adni_processed_fsl", recursive = FALSE)

lhipp <- matrix(ncol = 6)
rhipp <- matrix(ncol = 6)

lhipp_vol <- vector()
rhipp_vol <- vector()

img_id <- vector()
img_date <- vector()

lhipp_surface <- matrix(ncol = 6)
rhipp_surface <- matrix(ncol = 6)

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
    img_id <- c(img_id, patnos[dir_idx])
    img_date <- c(img_date, scan_dates[scan_idx])
    print(scan_dates[scan_idx])
    if (file.exists(paste0(scan_dirs[scan_idx], "-L_Hipp_first.nii.gz"))) {
      temp_lhipp <- readnii(paste0(scan_dirs[scan_idx], "-L_Hipp_first.nii.gz"))
      temp_lhipp_vol <- fslstats(
        paste0(scan_dirs[scan_idx], "_all_fast_firstseg.nii.gz"),
        opts = "-l 16.5 -u 17.5 -V"
      ) %>%
        str_replace(
          pattern = "([^\\s]+)",
          replacement = ""
        ) %>%
        gsub(
          pattern = " ",
          replacement = ""
        ) %>%
        as.numeric()
      lhipp_vol <- c(lhipp_vol, temp_lhipp_vol)
      lhipp_idx <- which(temp_lhipp > 0, arr.ind = TRUE)
      surface <- rep(FALSE, nrow(lhipp_idx))
      for (dim in 1:3) {
        for (dim2 in 1:3) {
          if (dim != dim2) {
            dim3 <- seq(1, 3, 1)[!(1:3 %in% c(dim, dim2))]
            unique_vals <- unique(lhipp_idx[, c(dim, dim2)])
            for (i in 1:nrow(unique_vals)) {
              vals <- lhipp_idx[lhipp_idx[, dim] == unique_vals[i, 1] & lhipp_idx[, dim2] == unique_vals[i, 2], dim3]
              min_vals <- vector()
              max_vals <- vector()
              for (v in vals) {
                if (!((v-1) %in% vals)) {
                  min_vals <- c(min_vals, v)
                } else if (!((v+1) %in% vals)) {
                  max_vals <- c(max_vals, v)
                }
              }
              surface_vals <- c(min_vals, max_vals)
              for (v in surface_vals) {
                surface[which(lhipp_idx[, dim] == unique_vals[i, 1] & lhipp_idx[, dim2] == unique_vals[i, 2] & lhipp_idx[, dim3] == v, arr.ind = TRUE)] <- TRUE
              }
            }
          } else {
            next
          }
        }
      }
      lhipp_temp <- cbind(
        patnos[dir_idx],
        scan_dates[scan_idx],
        lhipp_idx %*% diag(pixdim(temp_lhipp)[2:4]),
        temp_lhipp[lhipp_idx]
      )
      lhipp_surface_temp <- cbind(
        patnos[dir_idx],
        scan_dates[scan_idx],
        lhipp_idx[surface, ] %*% diag(pixdim(temp_lhipp)[2:4]),
        temp_lhipp[lhipp_idx[surface, ]]
      )
      lhipp <- rbind(lhipp, lhipp_temp)
      lhipp_surface <- rbind(lhipp_surface, lhipp_surface_temp)
    }
    if (file.exists(paste0(scan_dirs[scan_idx], "-R_Hipp_first.nii.gz"))) {
      temp_rhipp <- readnii(paste0(scan_dirs[scan_idx], "-R_Hipp_first.nii.gz"))
      temp_rhipp_vol <- fslstats(
        paste0(scan_dirs[scan_idx], "_all_fast_firstseg.nii.gz"),
        opts = "-l 52.5 -u 53.5 -V"
      ) %>%
        str_replace(
          pattern = "([^\\s]+)",
          replacement = ""
        ) %>%
        gsub(
          pattern = " ",
          replacement = ""
        ) %>%
        as.numeric()
      rhipp_vol <- c(rhipp_vol, temp_rhipp_vol)
      rhipp_idx <- which(temp_rhipp > 0, arr.ind = TRUE)
      surface <- rep(FALSE, nrow(rhipp_idx))
      for (dim in 1:3) {
        for (dim2 in 1:3) {
          if (dim != dim2) {
            dim3 <- seq(1, 3, 1)[!(1:3 %in% c(dim, dim2))]
            unique_vals <- unique(rhipp_idx[, c(dim, dim2)])
            for (i in 1:nrow(unique_vals)) {
              vals <- rhipp_idx[rhipp_idx[, dim] == unique_vals[i, 1] & rhipp_idx[, dim2] == unique_vals[i, 2], dim3]
              min_vals <- vector()
              max_vals <- vector()
              for (v in vals) {
                if (!((v-1) %in% vals)) {
                  min_vals <- c(min_vals, v)
                } else if (!((v+1) %in% vals)) {
                  max_vals <- c(max_vals, v)
                }
              }
              surface_vals <- c(min_vals, max_vals)
              for (v in surface_vals) {
                surface[which(rhipp_idx[, dim] == unique_vals[i, 1] & rhipp_idx[, dim2] == unique_vals[i, 2] & rhipp_idx[, dim3] == v, arr.ind = TRUE)] <- TRUE
              }
            }
          } else {
            next
          }
        }
      }
      rhipp_temp <- cbind(
        patnos[dir_idx],
        scan_dates[scan_idx],
        rhipp_idx %*% diag(pixdim(temp_rhipp)[2:4]),
        temp_rhipp[rhipp_idx]
      )
      rhipp_surface_temp <- cbind(
        patnos[dir_idx],
        scan_dates[scan_idx],
        rhipp_idx[surface, ] %*% diag(pixdim(temp_rhipp)[2:4]),
        temp_rhipp[rhipp_idx[surface, ]]
      )
      rhipp <- rbind(rhipp, rhipp_temp)
      rhipp_surface <- rbind(rhipp_surface, rhipp_surface_temp)
    }
  }
}

lhipp <- lhipp[-1, ] %>%
  as_tibble() %>%
  rename(
    patno = V1,
    scan_date = V2,
    x = V3,
    y = V4,
    z = V5,
    intensity = V6
  ) %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y),
    z = as.numeric(z),
    intensity = as.numeric(intensity)
  )

lhipp_surface <- lhipp_surface[-1, ] %>%
  as_tibble() %>%
  rename(
    patno = V1,
    scan_date = V2,
    x = V3,
    y = V4,
    z = V5,
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
    x = V3,
    y = V4,
    z = V5,
    intensity = V6
  ) %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y),
    z = as.numeric(z),
    intensity = as.numeric(intensity)
  )

rhipp_surface <- rhipp_surface[-1, ] %>%
  as_tibble() %>%
  rename(
    patno = V1,
    scan_date = V2,
    x = V3,
    y = V4,
    z = V5,
    intensity = V6
  ) %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y),
    z = as.numeric(z),
    intensity = as.numeric(intensity)
  )

hipp_info <- data.frame(img_id, img_date, lhipp_vol, rhipp_vol)
names(hipp_info) <- c("patno", "date", "lhipp_vol", "rhipp_vol")

write_csv(hipp_info, "data/adni_fsl_hipp_info.csv")

write_csv(lhipp, "data/adni_fsl_lhipp.csv")
write_csv(rhipp, "data/adni_fsl_rhipp.csv")

write_csv(lhipp_surface, "data/adni_fsl_lhipp_surface.csv")
write_csv(rhipp_surface, "data/adni_fsl_rhipp_surface.csv")
