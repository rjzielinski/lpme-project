library(pracma)
library(lubridate)
library(tidyverse)
source("code/pme.R")
source("code/lpme.R")

lhipp_surface <- read_csv("data/adni_fsl_lhipp_surface_red.csv")

lhipp_surface <- lhipp_surface %>%
  mutate(
    scan_date = gsub("\\.0", "", scan_date),
    scan_date = ymd_hms(scan_date),
    scan_date = decimal_date(scan_date)
  )

lhipp_surface_centers <- lhipp_surface %>%
  group_by(
    patno,
    scan_date
  ) %>%
  summarize(
    mean_x = mean(x),
    mean_y = mean(y),
    mean_z = mean(z),
    max_x = max(abs(x)),
    max_y = max(abs(y)),
    max_z = max(abs(z))
  )

lhipp_bl <- lhipp_surface %>%
  group_by(patno) %>%
  arrange(scan_date) %>%
  summarize(
    time_bl = first(scan_date),
    max_time = max(scan_date)
  ) %>%
  mutate(duration = max_time - time_bl)

overlap <- 0.05

lhipp_surface <- lhipp_surface %>%
  full_join(lhipp_surface_centers, by = c("patno", "scan_date")) %>%
  full_join(lhipp_bl, by = "patno") %>%
  mutate(
    x = x - mean_x,
    y = y - mean_y,
    z = z - mean_z,
    range_x = max(x) - min(x),
    range_y = max(y) - min(y),
    range_z = max(z) - min(z)
  ) %>%
  mutate(
    partition1 = ifelse((x > (0 - (overlap * range_x))) & (y > (0 - (overlap * range_y))), TRUE, FALSE),
    partition2 = ifelse((x > (0 - (overlap * range_x))) & (y <= (0 + (overlap * range_y))), TRUE, FALSE),
    partition3 = ifelse((x <= (0 + (overlap * range_x))) & (y > (0 - (overlap * range_y))), TRUE, FALSE),
    partition4 = ifelse((x <= (0 + (overlap * range_x))) & (y <= (0 + (overlap * range_y))), TRUE, FALSE),
    partition_sum = partition1 + partition2 + partition3 + partition4,
    partition = case_when(
      partition1 == TRUE & partition_sum == 1 ~ 1,
      partition2 == TRUE & partition_sum == 1 ~ 2,
      partition3 == TRUE & partition_sum == 1 ~ 3,
      partition4 == TRUE & partition_sum == 1 ~ 4,
      partition_sum > 1 ~ 5
    ),
    time_from_bl = scan_date - time_bl
    # r = sqrt(x^2 + y^2 + z^2),
    # theta = atan(y / x),
    # phi = acos(z / r)
  ) %>%
  mutate(
    x = x / max_x,
    y = y / max_y,
    z = z / max_z,
    time_from_bl = time_from_bl / duration
  )

lhipp_surface_spherical <- lhipp_surface %>%
  dplyr::select(x, y, z) %>%
  as.matrix() %>%
  cart2sph() %>%
  as_tibble()

lhipp_surface <- bind_cols(
  lhipp_surface,
  lhipp_surface_spherical
)

lhipp_test <- lhipp_surface %>%
  filter(patno == "002_S_0413")

lhipp_test_pt1 <- lhipp_test %>%
  filter(partition1 == TRUE)

lhipp_test_pt2 <- lhipp_test %>%
  filter(partition2 == TRUE)

lhipp_test_pt3 <- lhipp_test %>%
  filter(partition3 == TRUE)

lhipp_test_pt4 <- lhipp_test %>%
  filter(partition4 == TRUE)

lhipp_test_mat <- lhipp_test %>%
  dplyr::select(
    time_from_bl,
    x,
    y,
    z,
    theta,
    phi,
    r
  ) %>%
  as.matrix()

lhipp_test_bl <- lhipp_test_mat[lhipp_test_mat[, 1] == 0, -1]

lhipp_pt1_mat <- lhipp_test_pt1 %>%
  dplyr::select(
    time_from_bl,
    x,
    y,
    z
  ) %>%
  as.matrix()

lhipp_pt2_mat <- lhipp_test_pt2 %>%
  dplyr::select(
    time_from_bl,
    x,
    y,
    z
  ) %>%
  as.matrix()

lhipp_pt3_mat <- lhipp_test_pt3 %>%
  dplyr::select(
    time_from_bl,
    x,
    y,
    z
  ) %>%
  as.matrix()

lhipp_pt4_mat <- lhipp_test_pt4 %>%
  dplyr::select(
    time_from_bl,
    x,
    y,
    z
  ) %>%
  as.matrix()

test_pme_pt1 <- pme(lhipp_pt1_mat[lhipp_pt1_mat[, 1] == 0, -1], 2)

lpme_test_pt1 <- lpme(lhipp_pt1_mat, 2)
lpme_test_pt2 <- lpme(lhipp_pt2_mat, 2)
lpme_test_pt3 <- lpme(lhipp_pt3_mat, 2)
lpme_test_pt4 <- lpme(lhipp_pt4_mat, 2)

test_pme <- pme(lhipp_test_bl, d = 2, print_plots = TRUE)
test_lpme <- lpme(lhipp_test_mat, d = 2, tuning.para.seq = exp(seq(-20, 5, 0.25)))


### Next step: glue estimated manifolds together to estimate full surface
