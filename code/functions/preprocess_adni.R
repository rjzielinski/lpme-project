preprocess_adni <- function(
  adni_surface,
  adni_info,
  min_duration = 2,
  min_visits = 3,
  scan_exclude_threshold = 0.0025
) {
  require(dplyr)
  require(lubridate)
  require(tidyr)

  adni_surface <- adni_surface |>
    left_join(adni_info, by = join_by(scan_id == image_id))

  adni_surface <- adni_surface |>
    mutate(date = decimal_date(date))

  adni_n <- adni_surface |>
    group_by(subid, scan_id) |>
    tally() |>
    ungroup()
  adni_threshold <- quantile(adni_n$n, scan_exclude_threshold)

  exclude_scans <- adni_n |>
    filter(n < adni_threshold) |>
    select(scan_id) |>
    unlist()

  adni_surface <- adni_surface |>
    filter(!(scan_id %in% exclude_scans))

  adni_centers <- adni_surface |>
    group_by(subid, date, scan_id) |>
    summarise(
      mean_x = mean(x),
      mean_y = mean(y),
      mean_z = mean(z),
      max_x = max(abs(x - mean_x)),
      max_y = max(abs(y - mean_y)),
      max_z = max(abs(z - mean_z))
    ) |>
    ungroup()

  adni_bl <- adni_surface |>
    group_by(subid) |>
    arrange(date) |>
    summarize(
      date_bl = first(date),
      date_max = max(date),
      n_visits = length(unique(date))
    ) |>
    mutate(duration = date_max - date_bl)

  adni_surface <- adni_surface |>
    left_join(adni_centers, by = c("subid", "date", "scan_id")) |>
    left_join(adni_bl, by = "subid") |>
    mutate(
      x = (x - mean_x) / max_x,
      y = (y - mean_y) / max_y,
      z = (z - mean_z) / max_z,
      time_from_bl = (date - date_bl) / duration
    ) |>
    filter(duration > min_duration, n_visits >= min_visits)

  adni_surface_spherical <- adni_surface |>
    select(x, y, z) |>
    as.matrix() |>
    cart2sph() |>
    as_tibble()

  adni_surface <- bind_cols(adni_surface, adni_surface_spherical)

  list(surface = adni_surface, centers = adni_centers)
}
