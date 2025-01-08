# function for creating cross-section plots
create_cross_section_matrix <- function(mats, slice_idx, slice_widths, nrow = 1, lhipp) {
  time_vals <- unlist(mats[[1]][, 1]) %>% 
    unique()
  slices <- list()
  for (time_idx in 1:length(time_vals)) {
    slice_list <- list()
    for (idx in 1:length(mats)) {
      mat <- mats[[idx]]
      temp_mat <- mat[almost.equal(mat[, 1], time_vals[time_idx]), ] %>% 
        as.matrix()
      n_pixels <- length(unique(temp_mat[, slice_idx]))
      slice_mean <- mean(temp_mat[, slice_idx])
      # in the PME- and LPME-generated value matrices, voxel values are not all
      # consistent as in the data
      # to approximate a cross-section through the middle of the structure, we
      # consider all voxel values within a certain range of the mean voxel value.
      slice_list[[idx]] <- temp_mat[(temp_mat[, slice_idx] >= slice_mean - slice_widths[idx]) &
                                      (temp_mat[, slice_idx] <= slice_mean + slice_widths[idx]), ] %>%
        cbind(idx)
    }
    slices[[time_idx]] <- reduce(slice_list, rbind)
  }
  mat_sliced <- reduce(slices, rbind)

  slice_df <- data.frame(mat_sliced)
  names(slice_df) <- c("time", "x", "y", "z", "theta", "phi", "r", "Source")
  slice_df <- slice_df %>%
    mutate(
      # time = time * unique(lhipp$duration),
      time = round(time, 2),
      time = as.factor(time),
      Source = ifelse(
        Source == 1,
        "Data",
        ifelse(
          Source == 2,
          "LPME",
          "PME"
        )
      ),
      # code to change opacity by data source if needed
      opacity = case_when(
        Source == "Data" ~ 1,
        Source == "LPME" ~ 1,
        Source == "PME" ~ 1
      )
    )

  slice_df <- arrange(slice_df, Source)
  slice_df <- slice_df %>%
    mutate(
      Source = reorder(
        Source,
        X = c(
          rep(1, nrow(mat_sliced[mat_sliced[, 8] == 1, ])),
          rep(3, nrow(mat_sliced[mat_sliced[, 8] == 2, ])),
          rep(2, nrow(mat_sliced[mat_sliced[, 8] == 3, ]))
        )
      )
    )

  indices <- c("x", "y", "z")[2:4 != slice_idx]
  lim_val <- max(max(abs(slice_df[[indices[1]]])), max(abs(slice_df[[indices[2]]])))
  plt <- ggplot(arrange(slice_df, desc(as.character(Source)))) +
    geom_point(
      aes(
        x = .data[[indices[1]]],
        y = .data[[indices[2]]],
        color = Source,
        shape = Source
      ),
      size = 1
    ) +
    xlim(-lim_val, lim_val) +
    ylim(-lim_val, lim_val) +
    facet_wrap(vars(time), nrow = nrow) +
    theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0))

  plt
}