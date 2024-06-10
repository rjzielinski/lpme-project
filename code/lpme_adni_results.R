library(tidyverse)
library(brolgar)
library(RColorBrewer)

colors <- brewer.pal(3, "Set1")

roster <- read_csv("data/adni_standard/ROSTER.csv")
roster <- roster %>% 
  select(
    -USERDATE,
    -USERDATE2,
    -update_stamp
  )
adni_status <- read_csv("data/adni_standard/ADNIMERGE.csv")

# adni_status <- full_join(
#   adni_status,
#   roster,
#   by = c("Phase", "ID", "RID", "SITEID")
# )

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
  filter(max_date - date_bl > 2) %>%
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

lhipp_data_vol_pred <- vector()
lhipp_lpme_vol_pred <- vector()
lhipp_pme_vol_pred <- vector()
rhipp_data_vol_pred <- vector()
rhipp_lpme_vol_pred <- vector()
rhipp_pme_vol_pred <- vector()

for (patno_val in unique(hipp_info$patno)) {
  patno_data <- hipp_info %>%
    filter(patno == patno_val)
  lhipp_data_lm <- lm(lhipp_data_vol2 ~ date, data = patno_data)
  lhipp_data_pred <- predict(lhipp_data_lm, newdata = patno_data)
  lhipp_data_vol_pred <- c(lhipp_data_vol_pred, lhipp_data_pred)

  lhipp_lpme_lm <- lm(lhipp_vol_lpme2 ~ date, data = patno_data)
  lhipp_lpme_pred <- predict(lhipp_lpme_lm, newdata = patno_data)
  lhipp_lpme_vol_pred <- c(lhipp_lpme_vol_pred, lhipp_lpme_pred)

  lhipp_pme_lm <- lm(lhipp_vol_pme2 ~ date, data = patno_data)
  lhipp_pme_pred <- predict(lhipp_pme_lm, newdata = patno_data)
  lhipp_pme_vol_pred <- c(lhipp_pme_vol_pred, lhipp_pme_pred)

  rhipp_data_lm <- lm(rhipp_data_vol2 ~ date, data = patno_data)
  rhipp_data_pred <- predict(rhipp_data_lm, newdata = patno_data)
  rhipp_data_vol_pred <- c(rhipp_data_vol_pred, rhipp_data_pred)

  rhipp_lpme_lm <- lm(rhipp_vol_lpme2 ~ date, data = patno_data)
  rhipp_lpme_pred <- predict(rhipp_lpme_lm, newdata = patno_data)
  rhipp_lpme_vol_pred <- c(rhipp_lpme_vol_pred, rhipp_lpme_pred)

  rhipp_pme_lm <- lm(rhipp_vol_pme2 ~ date, data = patno_data)
  rhipp_pme_pred <- predict(rhipp_pme_lm, newdata = patno_data)
  rhipp_pme_vol_pred <- c(rhipp_pme_vol_pred, rhipp_pme_pred)
}

hipp_lm_pred <- tibble(
  lhipp_data_vol_pred,
  lhipp_lpme_vol_pred,
  lhipp_pme_vol_pred,
  rhipp_data_vol_pred,
  rhipp_lpme_vol_pred,
  rhipp_pme_vol_pred
)

hipp_info <- bind_cols(hipp_info, hipp_lm_pred)

hipp_info_sd <- hipp_info %>%
  group_by(patno) %>%
  summarize(
    lhipp_data_sd = sd(lhipp_data_vol2),
    lhipp_lpme_sd = sd(lhipp_vol_lpme2),
    lhipp_pme_sd = sd(lhipp_vol_pme2),
    rhipp_data_sd = sd(rhipp_data_vol2),
    rhipp_lpme_sd = sd(rhipp_vol_lpme2),
    rhipp_pme_sd = sd(rhipp_vol_pme2),
    lhipp_data_adj_sd = sd(lhipp_data_vol2 - lhipp_data_vol_pred),
    lhipp_lpme_adj_sd = sd(lhipp_vol_lpme2 - lhipp_lpme_vol_pred),
    lhipp_pme_adj_sd = sd(lhipp_vol_pme2 - lhipp_pme_vol_pred),
    rhipp_data_adj_sd = sd(rhipp_data_vol2 - rhipp_data_vol_pred),
    rhipp_lpme_adj_sd = sd(rhipp_vol_lpme2 - rhipp_lpme_vol_pred),
    rhipp_pme_adj_sd = sd(rhipp_vol_pme2 - rhipp_pme_vol_pred)
  ) %>%
  ungroup()

hipp_sd_mean <- hipp_info_sd %>%
  summarize(
    lhipp_data_sd_mean = mean(lhipp_data_sd),
    lhipp_lpme_sd_mean = mean(lhipp_lpme_sd),
    lhipp_pme_sd_mean = mean(lhipp_pme_sd),
    rhipp_data_sd_mean = mean(rhipp_data_sd),
    rhipp_lpme_sd_mean = mean(rhipp_lpme_sd),
    rhipp_pme_sd_mean = mean(rhipp_pme_sd),
    lhipp_data_adj_sd_mean = mean(lhipp_data_adj_sd),
    lhipp_lpme_adj_sd_mean = mean(lhipp_lpme_adj_sd),
    lhipp_pme_adj_sd_mean = mean(lhipp_pme_adj_sd),
    rhipp_data_adj_sd_mean = mean(rhipp_data_adj_sd),
    rhipp_lpme_adj_sd_mean = mean(rhipp_lpme_adj_sd),
    rhipp_pme_adj_sd_mean = mean(rhipp_pme_adj_sd)
)

print(hipp_sd_mean, width = Inf)

hipp_sd_med <- hipp_info_sd %>%
  summarize(
    lhipp_data_sd_med = median(lhipp_data_sd),
    lhipp_lpme_sd_med = median(lhipp_lpme_sd),
    lhipp_pme_sd_med = median(lhipp_pme_sd),
    rhipp_data_sd_med = median(rhipp_data_sd),
    rhipp_lpme_sd_med = median(rhipp_lpme_sd),
    rhipp_pme_sd_med = median(rhipp_pme_sd),
    lhipp_data_adj_sd_med = median(lhipp_data_adj_sd),
    lhipp_lpme_adj_sd_med = median(lhipp_lpme_adj_sd),
    lhipp_pme_adj_sd_med = median(lhipp_pme_adj_sd),
    rhipp_data_adj_sd_med = median(rhipp_data_adj_sd),
    rhipp_lpme_adj_sd_med = median(rhipp_lpme_adj_sd),
    rhipp_pme_adj_sd_med = median(rhipp_pme_adj_sd)
  )

print(hipp_sd_med, width = Inf)

thal_info_ts %>%
  sample_n_keys(3) %>%
  ggplot(aes(x = time_from_bl, group = patno, color = patno)) +
  geom_line(aes(y = lthal_data_vol2), linetype = "solid") +
  geom_line(aes(y = lthal_vol_lpme2), linetype = "dashed") +
  geom_line(aes(y = lthal_vol_pme2), linetype = "dotted") +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Thalamus Volume")

lthal_data_vol_pred <- vector()
lthal_lpme_vol_pred <- vector()
lthal_pme_vol_pred <- vector()
rthal_data_vol_pred <- vector()
rthal_lpme_vol_pred <- vector()
rthal_pme_vol_pred <- vector()

for (patno_val in unique(thal_info$patno)) {
  patno_data <- thal_info %>%
    filter(patno == patno_val)
  lthal_data_lm <- lm(lthal_data_vol2 ~ date, data = patno_data)
  lthal_data_pred <- predict(lthal_data_lm, newdata = patno_data)
  lthal_data_vol_pred <- c(lthal_data_vol_pred, lthal_data_pred)

  lthal_lpme_lm <- lm(lthal_vol_lpme2 ~ date, data = patno_data)
  lthal_lpme_pred <- predict(lthal_lpme_lm, newdata = patno_data)
  lthal_lpme_vol_pred <- c(lthal_lpme_vol_pred, lthal_lpme_pred)

  lthal_pme_lm <- lm(lthal_vol_pme2 ~ date, data = patno_data)
  lthal_pme_pred <- predict(lthal_pme_lm, newdata = patno_data)
  lthal_pme_vol_pred <- c(lthal_pme_vol_pred, lthal_pme_pred)

  rthal_data_lm <- lm(rthal_data_vol2 ~ date, data = patno_data)
  rthal_data_pred <- predict(rthal_data_lm, newdata = patno_data)
  rthal_data_vol_pred <- c(rthal_data_vol_pred, rthal_data_pred)

  rthal_lpme_lm <- lm(rthal_vol_lpme2 ~ date, data = patno_data)
  rthal_lpme_pred <- predict(rthal_lpme_lm, newdata = patno_data)
  rthal_lpme_vol_pred <- c(rthal_lpme_vol_pred, rthal_lpme_pred)

  rthal_pme_lm <- lm(rthal_vol_pme2 ~ date, data = patno_data)
  rthal_pme_pred <- predict(rthal_pme_lm, newdata = patno_data)
  rthal_pme_vol_pred <- c(rthal_pme_vol_pred, rthal_pme_pred)
}

thal_lm_pred <- tibble(
  lthal_data_vol_pred,
  lthal_lpme_vol_pred,
  lthal_pme_vol_pred,
  rthal_data_vol_pred,
  rthal_lpme_vol_pred,
  rthal_pme_vol_pred
)

thal_info <- bind_cols(thal_info, thal_lm_pred)


thal_info_sd <- thal_info %>%
  group_by(patno) %>%
  summarize(
    lthal_data_sd = sd(lthal_data_vol2),
    lthal_lpme_sd = sd(lthal_vol_lpme2),
    lthal_pme_sd = sd(lthal_vol_pme2),
    rthal_data_sd = sd(rthal_data_vol2),
    rthal_lpme_sd = sd(rthal_vol_lpme2),
    rthal_pme_sd = sd(rthal_vol_pme2),
    lthal_data_adj_sd = sd(lthal_data_vol2 - lthal_data_vol_pred),
    lthal_lpme_adj_sd = sd(lthal_vol_lpme2 - lthal_lpme_vol_pred),
    lthal_pme_adj_sd = sd(lthal_vol_pme2 - lthal_pme_vol_pred),
    rthal_data_adj_sd = sd(rthal_data_vol2 - rthal_data_vol_pred),
    rthal_lpme_adj_sd = sd(rthal_vol_lpme2 - rthal_lpme_vol_pred),
    rthal_pme_adj_sd = sd(rthal_vol_pme2 - rthal_pme_vol_pred)
  ) %>%
  ungroup()

thal_sd_mean <- thal_info_sd %>%
  summarize(
    lthal_data_sd_mean = mean(lthal_data_sd),
    lthal_lpme_sd_mean = mean(lthal_lpme_sd),
    lthal_pme_sd_mean = mean(lthal_pme_sd),
    rthal_data_sd_mean = mean(rthal_data_sd),
    rthal_lpme_sd_mean = mean(rthal_lpme_sd),
    rthal_pme_sd_mean = mean(rthal_pme_sd),
    lthal_data_adj_sd_mean = mean(lthal_data_adj_sd),
    lthal_lpme_adj_sd_mean = mean(lthal_lpme_adj_sd),
    lthal_pme_adj_sd_mean = mean(lthal_pme_adj_sd),
    rthal_data_adj_sd_mean = mean(rthal_data_adj_sd),
    rthal_lpme_adj_sd_mean = mean(rthal_lpme_adj_sd),
    rthal_pme_adj_sd_mean = mean(rthal_pme_adj_sd)
  )

print(thal_sd_mean, width = Inf)

thal_sd_med <- thal_info_sd %>%
  summarize(
    lthal_data_sd_med = median(lthal_data_sd),
    lthal_lpme_sd_med = median(lthal_lpme_sd),
    lthal_pme_sd_med = median(lthal_pme_sd),
    rthal_data_sd_med = median(rthal_data_sd),
    rthal_lpme_sd_med = median(rthal_lpme_sd),
    rthal_pme_sd_med = median(rthal_pme_sd),
    lthal_data_adj_sd_med = median(lthal_data_adj_sd),
    lthal_lpme_adj_sd_med = median(lthal_lpme_adj_sd),
    lthal_pme_adj_sd_med = median(lthal_pme_adj_sd),
    rthal_data_adj_sd_med = median(rthal_data_adj_sd),
    rthal_lpme_adj_sd_med = median(rthal_lpme_adj_sd),
    rthal_pme_adj_sd_med = median(rthal_pme_adj_sd)
  )

print(thal_sd_med, width = Inf)

set.seed(8310)
sampled_patnos <- hipp_info_ts %>% 
  sample_n_keys(3) %>% 
  .$patno %>% 
  unique()
hipp_info_ts %>%
  filter(patno %in% sampled_patnos) %>%
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_point(aes(y = lhipp_data_vol2), shape = 16, color = colors[1]) +
  geom_line(aes(y = lhipp_data_vol2), color = colors[1]) +
  # geom_smooth(
  #   aes(y = lhipp_data_vol2),
  #   color = colors[1],
  #   method = "lm",
  #   se = FALSE,
  #   linetype = "dashed",
  #   linewidth = 1
  # ) +
  geom_point(aes(y = lhipp_vol_lpme2), shape = 15, color = colors[2]) +
  geom_line(aes(y = lhipp_vol_lpme2), color = colors[2]) +
  # geom_smooth(
  #   aes(y = lhipp_vol_lpme2),
  #   color = colors[2],
  #   method = "lm",
  #   se = FALSE,
  #   linetype = "dashed",
  #   linewidth = 1
  # ) +
  geom_point(aes(y = lhipp_vol_pme2), shape = 17, color = colors[3]) +
  geom_line(aes(y = lhipp_vol_pme2), color = colors[3]) +
  # geom_smooth(
  #   aes(y = lhipp_vol_pme2),
  #   color = colors[3],
  #   method = "lm",
  #   se = FALSE,
  #   linetype = "dashed",
  #   linewidth = 1
  # ) +
  facet_wrap(~patno) +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Hippocampus Volume")
ggsave("paper/figures/adni_plots/adni_lhipp_volume_comp.png", dpi = 1500)

thal_info_ts %>%
  filter(patno %in% sampled_patnos) %>% 
  ggplot(aes(x = time_from_bl, group = patno)) +
  geom_point(aes(y = lthal_data_vol2), shape = 16, color = colors[1]) +
  geom_line(aes(y = lthal_data_vol2), linetype = "solid", color = colors[1]) +
  # geom_smooth(
  #   aes(y = lthal_data_vol2),
  #   color = colors[1],
  #   method = "lm",
  #   se = FALSE,
  #   linetype = "dashed",
  #   linewidth = 1
  # ) +
  geom_point(aes(y = lthal_vol_lpme2), shape = 15, color = colors[2]) +
  geom_line(aes(y = lthal_vol_lpme2), linetype = "dashed", color = colors[2]) +
  # geom_smooth(
  #   aes(y = lthal_vol_lpme2),
  #   color = colors[2],
  #   method = "lm",
  #   se = FALSE,
  #   linetype = "dashed",
  #   linewidth = 1
  # ) +
  geom_point(aes(y = lthal_vol_pme2), shape = 17, color = colors[3]) +
  geom_line(aes(y = lthal_vol_pme2), linetype = "dotted", color = colors[3]) +
  # geom_smooth(
  #   aes(y = lthal_vol_pme2),
  #   color = colors[3],
  #   method = "lm",
  #   se = FALSE,
  #   linetype = "dashed",
  #   linewidth = 1
  # ) +
  facet_wrap(~patno) +
  xlab("Time from Baseline Visit (Years)") +
  ylab("Estimated Left Thalamus Volume") + 
  guides(
    color = guide_legend(title = "Participant ID"),
    shape = guide_legend(title = "Estimate Source")
  )
ggsave("paper/figures/adni_plots/adni_lthal_volume_comp.png", dpi = 1500)

adni_status %>% 
  filter(
    PTID %in% eligible_patnos,
    VISCODE == "bl"
  ) %>% 
  group_by(DX_bl) %>% 
  tally()
