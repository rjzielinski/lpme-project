library(oro.dicom)
library(oro.nifti)
library(neurobase)
library(fslr)
library(extrantsr)
library(stringr)

sub_dirs <- list.dirs("data/ADNI", recursive = FALSE)

nii_scans <- list()

processed_scans <- list()

hipp_mask <- readNIfTI(
  paste0(fsl_dir(), "/data/standard/MNI152_T1_1mm_Hipp_mask_dil8.nii.gz")
)

for (dir_idx in 1:length(sub_dirs)) {
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

test_dicom <- readDICOM(
  list.files(sub_dirs[1], recursive = TRUE, full.names = TRUE)[1]
)
