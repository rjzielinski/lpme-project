library(oro.dicom)
library(oro.nifti)
library(neurobase)
library(fslr)
# library(extrantsr)
library(stringr)
library(foreach)
library(tidyverse)
library(ANTsRCore)

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
# foreach(dir_idx = 1:8) %dopar% {
  img_dirs <- list.files(sub_dirs[dir_idx], recursive = TRUE, full.names = TRUE) %>%
    gsub(pattern = "([^/]+$)", replacement = "") %>%
    unique()

  nii_scans[[dir_idx]] <- list()

  for (img_idx in 1:length(img_dirs)) {
    proc_dir <- paste0("data/adni_processed_fsl", gsub(pattern = "data/adni", replacement = "", img_dirs[img_idx]))
    if (file.exists(paste0(proc_dir, "/_all_fast_origsegs.nii.gz"))) {
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
      # nii <- bias_correct(
      #   nii,
      #   correction = "N4",
      #   verbose = FALSE
      # ) # N4 bias correction using ANTs
      # nii_bc <- fast(nii, bias_correct = TRUE)

      # nii <- fast(nii, bias_correct = TRUE)
      # nii <- n4BiasFieldCorrection(
      #   nii,
      #   verbose = FALSE
      # )

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
        verbose = FALSE,
        opts = "-d"
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

lthal <- matrix(ncol = 6)
rthal <- matrix(ncol = 6)

lhipp_vol <- vector()
rhipp_vol <- vector()

lthal_vol <- vector()
rthal_vol <- vector()

img_id <- vector()
img_date <- vector()

lhipp_surface <- matrix(ncol = 6)
rhipp_surface <- matrix(ncol = 6)

lthal_surface <- matrix(ncol = 6)
rthal_surface <- matrix(ncol = 6)

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
    if (file.exists(paste0(scan_dirs[scan_idx], "_all_fast_origsegs.nii.gz"))) {
      img_id <- c(img_id, patnos[dir_idx])
      img_date <- c(img_date, scan_dates[scan_idx])
      temp_segs <- readnii(paste0(scan_dirs[scan_idx], "_all_fast_origsegs.nii.gz"))
      seg_vols <- fslstats(
        paste0(scan_dirs[scan_idx], "_all_fast_origsegs.nii.gz"),
        opts = "-V",
        ts = TRUE
      )
      temp_lhipp <- temp_segs[, , , 6]
      temp_lhipp_vol <- seg_vols[6] %>%
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
        lhipp_idx %*% diag(pixdim(as.nifti(temp_lhipp))[2:4]),
        temp_lhipp[lhipp_idx]
      )
      lhipp_surface_temp <- cbind(
        patnos[dir_idx],
        scan_dates[scan_idx],
        lhipp_idx[surface, ] %*% diag(pixdim(as.nifti(temp_lhipp))[2:4]),
        temp_lhipp[lhipp_idx[surface, ]]
      )
      lhipp <- rbind(lhipp, lhipp_temp)
      lhipp_surface <- rbind(lhipp_surface, lhipp_surface_temp)
      temp_rhipp <- temp_segs[, , , 13]
      
      temp_rhipp_vol <- seg_vols[13] %>%
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
        rhipp_idx %*% diag(pixdim(as.nifti(temp_rhipp))[2:4]),
        temp_rhipp[rhipp_idx]
      )
      rhipp_surface_temp <- cbind(
        patnos[dir_idx],
        scan_dates[scan_idx],
        rhipp_idx[surface, ] %*% diag(pixdim(as.nifti(temp_rhipp))[2:4]),
        temp_rhipp[rhipp_idx[surface, ]]
      )
      rhipp <- rbind(rhipp, rhipp_temp)
      rhipp_surface <- rbind(rhipp_surface, rhipp_surface_temp)

      temp_lthal <- temp_segs[, , , 1]
      temp_lthal_vol <- seg_vols[1] %>%
        str_replace(
          pattern = "([^\\s]+)",
          replacement = ""
        ) %>%
        gsub(
          pattern = " ",
          replacement = ""
        ) %>%
        as.numeric()
      lthal_vol <- c(lthal_vol, temp_lthal_vol)
      lthal_idx <- which(temp_lthal > 0, arr.ind = TRUE)
      surface <- rep(FALSE, nrow(lthal_idx))
      for (dim in 1:3) {
        for (dim2 in 1:3) {
          if (dim != dim2) {
            dim3 <- seq(1, 3, 1)[!(1:3 %in% c(dim, dim2))]
            unique_vals <- unique(lthal_idx[, c(dim, dim2)])
            for (i in 1:nrow(unique_vals)) {
              vals <- lthal_idx[lthal_idx[, dim] == unique_vals[i, 1] & lthal_idx[, dim2] == unique_vals[i, 2], dim3]
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
                surface[which(lthal_idx[, dim] == unique_vals[i, 1] & lthal_idx[, dim2] == unique_vals[i, 2] & lthal_idx[, dim3] == v, arr.ind = TRUE)] <- TRUE
              }
            }
          } else {
            next
          }
        }
      }
      lthal_temp <- cbind(
        patnos[dir_idx],
        scan_dates[scan_idx],
        lthal_idx %*% diag(pixdim(as.nifti(temp_lthal))[2:4]),
        temp_lthal[lthal_idx]
      )
      lthal_surface_temp <- cbind(
        patnos[dir_idx],
        scan_dates[scan_idx],
        lthal_idx[surface, ] %*% diag(pixdim(as.nifti(temp_lthal))[2:4]),
        temp_lthal[lthal_idx[surface, ]]
      )
      lthal <- rbind(lthal, lthal_temp)
      lthal_surface <- rbind(lthal_surface, lthal_surface_temp)
      temp_rthal <- temp_segs[, , , 9]
      temp_rthal_vol <- seg_vols[9] %>%
        str_replace(
          pattern = "([^\\s]+)",
          replacement = ""
        ) %>%
        gsub(
          pattern = " ",
          replacement = ""
        ) %>%
        as.numeric()
      rthal_vol <- c(rthal_vol, temp_rthal_vol)
      rthal_idx <- which(temp_rthal > 0, arr.ind = TRUE)
      surface <- rep(FALSE, nrow(rthal_idx))
      for (dim in 1:3) {
        for (dim2 in 1:3) {
          if (dim != dim2) {
            dim3 <- seq(1, 3, 1)[!(1:3 %in% c(dim, dim2))]
            unique_vals <- unique(rthal_idx[, c(dim, dim2)])
            for (i in 1:nrow(unique_vals)) {
              vals <- rthal_idx[rthal_idx[, dim] == unique_vals[i, 1] & rthal_idx[, dim2] == unique_vals[i, 2], dim3]
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
                surface[which(rthal_idx[, dim] == unique_vals[i, 1] & rthal_idx[, dim2] == unique_vals[i, 2] & rthal_idx[, dim3] == v, arr.ind = TRUE)] <- TRUE
              }
            }
          } else {
            next
          }
        }
      }
      rthal_temp <- cbind(
        patnos[dir_idx],
        scan_dates[scan_idx],
        rthal_idx %*% diag(pixdim(as.nifti(temp_rthal))[2:4]),
        temp_rthal[rthal_idx]
      )
      rthal_surface_temp <- cbind(
        patnos[dir_idx],
        scan_dates[scan_idx],
        rthal_idx[surface, ] %*% diag(pixdim(as.nifti(temp_rthal))[2:4]),
        temp_rthal[rthal_idx[surface, ]]
      )
      rthal <- rbind(rthal, rthal_temp)
      rthal_surface <- rbind(rthal_surface, rthal_surface_temp)
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

lthal <- lthal[-1, ] %>%
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

lthal_surface <- lthal_surface[-1, ] %>%
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

rthal <- rthal[-1, ] %>%
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

rthal_surface <- rthal_surface[-1, ] %>%
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

thal_info <- data.frame(img_id, img_date, lthal_vol, rthal_vol)
names(thal_info) <- c("patno", "date", "lthal_vol", "rthal_vol")

write_csv(hipp_info, "data/adni_fsl_hipp_info.csv")
write_csv(thal_info, "data/adni_fsl_thal_info.csv")

write_csv(lhipp, "data/adni_fsl_lhipp.csv")
write_csv(rhipp, "data/adni_fsl_rhipp.csv")

write_csv(lthal, "data/adni_fsl_lthal.csv")
write_csv(rthal, "data/adni_fsl_rthal.csv")

write_csv(lhipp_surface, "data/adni_fsl_lhipp_surface.csv")
write_csv(rhipp_surface, "data/adni_fsl_rhipp_surface.csv")

write_csv(lthal_surface, "data/adni_fsl_lthal_surface.csv")
write_csv(rthal_surface, "data/adni_fsl_rthal_surface.csv")
