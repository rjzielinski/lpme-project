library(tidyverse)
library(brolgar)
library(RColorBrewer)

colors <- brewer.pal(3, "Set1")

sub_dirs <- list.dirs("results", recursive = FALSE)

hipp_info <- data.frame(
  patno = character(),
  date = numeric(),
  lhipp_data_vol2 = numeric(),
  lhipp_vol_lpme1 = numeric(),
  lhipp_vol_lpme2 = numeric(),
  lhipp_vol_pme1 = numeric(),
  lhipp_vol_pme2 = numeric(),
  rhipp_data_vol2 = numeric(),
  rhipp_vol_lpme1 = numeric(),
  rhipp_vol_lpme2 = numeric(),
  rhipp_vol_pme1 = numeric(),
  rhipp_vol_pme2 = numeric()
)

thal_info <- data.frame(
  patno = character(),
  date = numeric(),
  lthal_data_vol2 = numeric(),
  lthal_vol_lpme1 = numeric(),
  lthal_vol_lpme2 = numeric(),
  lthal_vol_pme1 = numeric(),
  lthal_vol_pme2 = numeric(),
  rthal_data_vol2 = numeric(),
  rthal_vol_lpme1 = numeric(),
  rthal_vol_lpme2 = numeric(),
  rthal_vol_pme1 = numeric(),
  rthal_vol_pme2 = numeric()
)

for (dir in sub_dirs) {
  file_list <- list.files(dir, full.names = FALSE)
  if ("hipp_info.csv" %in% file_list) {
    temp_hipp_info <- read_csv(paste0(dir, "/hipp_info.csv"))[, -1]
    hipp_info <- bind_rows(hipp_info, temp_hipp_info)
  }
  if ("thal_info.csv" %in% file_list) {
    temp_thal_info <- read_csv(paste0(dir, "/thal_info.csv"))[, -1]
    thal_info <- bind_rows(thal_info, temp_thal_info)
  }
}

hipp_bl <- hipp_info %>%
  group_by(patno) %>%
  summarize(
    date_bl = min(date),
    max_date = max(date)
  )
thal_bl <- thal_info %>%
  group_by(patno) %>%
  summarize(
    date_bl = min(date),
    max_date = max(date)
  )

eligible_patnos <- hipp_bl %>%
  filter(max_date - date_bl > 3) %>%
  select(patno) %>%
  unlist() %>%
  unique()

hipp_info <- full_join(hipp_info, hipp_bl, by = "patno") %>%
  mutate(time_from_bl = date - date_bl) %>%
  filter(patno %in% eligible_patnos)

thal_info <- full_join(thal_info, thal_bl, by = "patno") %>%
  mutate(time_from_bl = date - date_bl) %>%
  filter(patno %in% eligible_patnos)

hipp_summary <- hipp_info %>%
  group_by(patno) %>%
  summarize(
    lhipp_data_vol2_sd = sd(lhipp_data_vol2),
    lhipp_vol_lpme1_sd = sd(lhipp_vol_lpme1),
    lhipp_vol_lpme2_sd = sd(lhipp_vol_lpme2),
    lhipp_vol_pme1_sd = sd(lhipp_vol_pme1),
    lhipp_vol_pme2_sd = sd(lhipp_vol_pme2),
    rhipp_data_vol2_sd = sd(rhipp_data_vol2),
    rhipp_vol_lpme1_sd = sd(rhipp_vol_lpme1),
    rhipp_vol_lpme2_sd = sd(rhipp_vol_lpme2),
    rhipp_vol_pme1_sd = sd(rhipp_vol_pme1),
    rhipp_vol_pme2_sd = sd(rhipp_vol_pme2)
  )

hipp_long <- hipp_info %>%
  pivot_longer(
    lhipp_data_vol2:lhipp_vol_pme2,
    names_to = "lhipp_method",
    values_to = "lhipp_volume"
  ) %>%
  pivot_longer(
    rhipp_data_vol2:rhipp_vol_pme2,
    names_to = "rhipp_method",
    values_to = "rhipp_volume"
  ) %>%
  mutate(
    lhipp_method = gsub("lhipp_", "", lhipp_method),
    rhipp_method = gsub("rhipp_", "", rhipp_method)
  ) %>%
  filter(lhipp_method == rhipp_method) %>%
  select(-rhipp_method) %>%
  separate(lhipp_method, c("source", "method"), sep = "_") %>%
  mutate(
    source_mod = ifelse(source == "data", source, method),
    method_mod = ifelse(source == "vol", source, method)
  ) %>%
  select(-source, -method) %>%
  rename(source = source_mod, method = method_mod)

hipp_long$method <- ifelse(grepl("pme1", hipp_long$source), "vol1", "vol2")
hipp_long$source <- gsub("pme1", "pme", hipp_long$source)
hipp_long$source <- gsub("pme2", "pme", hipp_long$source)

thal_long <- thal_info %>%
  pivot_longer(
    lthal_data_vol2:lthal_vol_pme2,
    names_to = "lthal_method",
    values_to = "lthal_volume"
  ) %>%
  pivot_longer(
    rthal_data_vol2:rthal_vol_pme2,
    names_to = "rthal_method",
    values_to = "rthal_volume"
  ) %>%
  mutate(
    lthal_method = gsub("lthal_", "", lthal_method),
    rthal_method = gsub("rthal_", "", rthal_method)
  ) %>%
  filter(lthal_method == rthal_method) %>%
  select(-rthal_method) %>%
  separate(lthal_method, c("source", "method"), sep = "_") %>%
  mutate(
    source_mod = ifelse(source == "data", source, method),
    method_mod = ifelse(source == "vol", source, method)
  ) %>%
  select(-source, -method) %>%
  rename(source = source_mod, method = method_mod)

thal_long$method <- ifelse(grepl("pme1", thal_long$source), "vol1", "vol2")
thal_long$source <- gsub("pme1", "pme", thal_long$source)
thal_long$source <- gsub("pme2", "pme", thal_long$source)

thal_info_ts <- as_tsibble(
  thal_info,
  key = patno,
  index = time_from_bl,
  regular = FALSE
)

hipp_info_ts <- as_tsibble(
  hipp_info,
  key = patno,
  index = time_from_bl,
  regular = FALSE
)


hipp_info_ts %>%
  sample_n_keys(3) %>%
  ggplot(aes(x = time_from_bl, group = patno, color = patno)) +
  geom_line(aes(y = lhipp_data_vol2), linetype = "solid") +
  geom_line(aes(y = lhipp_vol_lpme2), linetype = "dashed") +
  geom_line(aes(y = lhipp_vol_pme2), linetype = "dotted") +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Hippocampus Volume")

hipp_info_sd <- hipp_info %>%
  group_by(patno) %>%
  summarize(
    data_sd = sd(lhipp_data_vol2),
    lpme_sd = sd(lhipp_vol_lpme2),
    pme_sd = sd(lhipp_vol_pme2),
    data_iqr = quantile(lhipp_data_vol2, 0.75) - quantile(lhipp_data_vol2, 0.25),
    lpme_iqr = quantile(lhipp_vol_lpme2, 0.75) - quantile(lhipp_vol_lpme2, 0.25),
    pme_iqr = quantile(lhipp_vol_pme2, 0.75) - quantile(lhipp_vol_pme2, 0.25)
  ) %>%
  ungroup()

hipp_info_sd %>%
  summarize(
    data_sd_mean = mean(data_sd),
    lpme_sd_mean = mean(lpme_sd),
    pme_sd_mean = mean(pme_sd),
    data_iqr_mean = mean(data_iqr),
    lpme_iqr_mean = mean(lpme_iqr),
    pme_iqr_mean = mean(pme_iqr)
  )

hipp_info_sd %>%
  summarize(
    data_sd_med = median(data_sd),
    lpme_sd_med = median(lpme_sd),
    pme_sd_med = median(pme_sd),
    data_iqr_med = median(data_iqr),
    lpme_iqr_med = median(lpme_iqr),
    pme_iqr_med = median(pme_iqr)
  )

thal_info_ts %>%
  sample_n_keys(3) %>%
  ggplot(aes(x = time_from_bl, group = patno, color = patno)) +
  geom_line(aes(y = lthal_data_vol2), linetype = "solid") +
  geom_line(aes(y = lthal_vol_lpme2), linetype = "dashed") +
  geom_line(aes(y = lthal_vol_pme2), linetype = "dotted") +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Thalamus Volume")

thal_info_sd <- thal_info %>%
  group_by(patno) %>%
  summarize(
    data_sd = sd(lthal_data_vol2),
    lpme_sd = sd(lthal_vol_lpme2),
    pme_sd = sd(lthal_vol_pme2),
    data_iqr = quantile(lthal_data_vol2, 0.75) - quantile(lthal_data_vol2, 0.25),
    lpme_iqr = quantile(lthal_vol_lpme2, 0.75) - quantile(lthal_vol_lpme2, 0.25),
    pme_iqr = quantile(lthal_vol_pme2, 0.75) - quantile(lthal_vol_pme2, 0.25)
  ) %>%
  ungroup()

thal_info_sd %>%
  summarize(
    data_sd_mean = mean(data_sd),
    lpme_sd_mean = mean(lpme_sd),
    pme_sd_mean = mean(pme_sd),
    data_iqr_mean = mean(data_iqr),
    lpme_iqr_mean = mean(lpme_iqr),
    pme_iqr_mean = mean(pme_iqr)
  )

thal_info_sd %>%
  summarize(
    data_sd_med = median(data_sd),
    lpme_sd_med = median(lpme_sd),
    pme_sd_med = median(pme_sd),
    data_iqr_med = median(data_iqr),
    lpme_iqr_med = median(lpme_iqr),
    pme_iqr_med = median(pme_iqr)
  )


set.seed(6819)
hipp_info_ts %>%
  sample_n_keys(6) %>%
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_line(aes(y = lhipp_data_vol2), color = colors[1]) +
  geom_line(aes(y = lhipp_vol_lpme2), color = colors[2]) +
  geom_line(aes(y = lhipp_vol_pme2), color = colors[3]) +
  facet_wrap(~ patno) +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Hippocampus Volume")
ggsave("paper/figures/adni_plots/adni_lhipp_volume_comp.png")

thal_info_ts %>%
  sample_n_keys(6) %>%
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_line(aes(y = lthal_data_vol2), color = colors[1]) +
  geom_line(aes(y = lthal_vol_lpme2), color = colors[2]) +
  geom_line(aes(y = lthal_vol_pme2), color = colors[3]) +
  facet_wrap(~ patno) +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Thalamus Volume")
ggsave("paper/figures/adni_plots/adni_lthal_volume_comp.png")


thal_long %>%
  filter(method == "vol2") %>%
  sample_n_keys(5) %>%
  ggplot(aes(x = date, y = lthal_volume, group = patno, color = source)) +
  geom_point() +
  geom_line()
