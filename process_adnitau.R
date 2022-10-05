library(oro.dicom)
library(oro.nifti)
library(tidyverse)
library(tools)
library(svMisc)

adni_tau_path <- "data/ADNITau/"

id_dirs <- list.files(adni_tau_path, full.names = TRUE)

for (dir_idx in 1:length(id_dirs)) {
  id_dir <- id_dirs[dir_idx]
  subdirs <- list.dirs(id_dir)
  num_nested <- grep("//", subdirs)
  deepest <- subdirs[num_nested == max(num_nested)]
  for (subdir in deepest) {
    file_names <- list.files(subdir, full.names = TRUE)
    if (sum(file_ext(file_names) == "dcm") == length(file_names)) {
      all_slices <- readDICOM(subdir)
      nii <- tryCatch(
        dicom2nifti(all_slices),
        error = function(e) {
          print(paste0("ERROR in directory: ", id_dir))
          print("")
        }
      )
      if (class(nii) != "character") {
        fname <- subdir %>% 
          gsub(
            pattern = "/",
            replacement = "_"
          ) %>% 
          gsub(
            pattern = "data_",
            replacement = ""
          ) %>% 
          gsub(
            pattern = "ADNITau",
            replacement = "ADNI"
          ) %>% 
          gsub(
            pattern = "__",
            replacement = "_"
          ) %>% 
          gsub(
            pattern = "\\.",
            replacement = "_"
          )
        
        writeNIfTI(
          nim = nii, 
          filename = paste(subdir, fname, sep = "/")
        )
        
        file.remove(file_names)
      }
    }
  }
  progress(dir_idx, max.value = length(id_dirs))
}

