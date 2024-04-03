library(tidyverse)
library(reticulate)
library(pme)
library(stringr)

# files_case1 <- list.files("simulations/case1", full.names = TRUE)
# files_case2 <- list.files("simulations/case2", full.names = TRUE)
# files_case3 <- list.files("simulations/case3", full.names = TRUE)
# files_case5 <- list.files("simulations/case5", full.names = TRUE)
# files_case6 <- list.files("simulations/case6", full.names = TRUE)
# files_case7 <- list.files("simulations/case7", full.names = TRUE)
# files_case8 <- list.files("simulations/case8", full.names = TRUE)
# files_case9 <- list.files("simulations/case9", full.names = TRUE)
# files_case10 <- list.files("simulations/case10", full.names = TRUE)
#
# files_full <- c(
#   files_case1,
#   files_case2,
#   files_case3,
#   files_case5,
#   files_case6,
#   files_case7,
#   files_case8,
#   files_case9,
#   files_case10
# )

files_full <- list.files("simulations", recursive = TRUE, full.names = TRUE)

# files_full <- files_case1

initialization_algorithm <- vector()
lpme_error <- vector()
pme_error <- vector()
# PME_error <- vector()
pme_ssd <- vector()
# PME_ssd <- vector()
principal_curve_ssd <- vector()
principal_curve_error <- vector()
data_error <- vector()
amp_noise <- vector()
period_noise <- vector()
n_visits <- vector()
max_time <- vector()
time_trend <- vector()
sim_case <- vector()
pme_ssd_raw <- vector()
principal_curve_ssd_raw <- vector()
for (i in 1:length(files_full)) {
  error_list <- readRDS(files_full[i])
  initialization_algorithm[i] <- error_list$initialization_algorithm
  lpme_error[i] <- error_list$lpme_error
  pme_error[i] <- error_list$pme_error
  # PME_error[i] <- error_list$PME_error
  principal_curve_error[i] <- error_list$principal_curve_error

  # PME_ssd <- c(PME_ssd, error_list$PME_ssd)
  principal_curve_ssd <- c(principal_curve_ssd, error_list$principal_curve_ssd)
  data_error[i] <- error_list$data_error
  amp_noise[i] <- error_list$amp_noise
  period_noise[i] <- error_list$period_noise
  n_visits[i] <- length(error_list$times)
  max_time[i] <- max(error_list$times)
  if (str_detect(files_full[i], "constant")) {
    time_trend[i] <- "constant"
  } else if (str_detect(files_full[i], "linear")) {
    time_trend[i] <- "linear"
  } else if (str_detect(files_full[i], "quadratic")) {
    time_trend[i] <- "quadratic"
  } else if (str_detect(files_full[i], "sinusoidal")) {
    time_trend[i] <- "sinusoidal"
  }
  case <- str_sub(files_full[i], 13, 18) %>%
    gsub(pattern = "/", replacement = "") %>%
    gsub(pattern = "case", replacement = "") %>%
    as.numeric()
  sim_case[i] <- case

  pme_ssd_vals <- map(
    error_list$pme_results,
    ~ min(.x$MSD) * (nrow(error_list$df) / n_visits[i])
  ) %>%
    reduce(c)
  pme_ssd_raw <- c(pme_ssd_raw, pme_ssd_vals)
  pme_ssd[i] <- mean(pme_ssd_vals)
  if (case < 7) {
    principal_curve_ssd_vals <-map(
      error_list$principal_curve_results,
      ~ .x$dist
    ) %>%
      reduce(c)
    principal_curve_ssd[i] <- mean(principal_curve_ssd_vals)
    principal_curve_ssd_raw <- c(principal_curve_ssd_raw, principal_curve_ssd_vals)
  } else {
    principal_curve_ssd_vals <- map(
      error_list$principal_curve_results,
      ~ .x$MSE * (nrow(error_list$df) / n_visits[i])
    ) %>%
      reduce(c)
    principal_curve_ssd[i] <- mean(principal_curve_ssd_vals)
    principal_curve_ssd_raw <- c(principal_curve_ssd_raw, principal_curve_ssd_vals)
  }
}

error_df <- data.frame(
  sim_case,
  initialization_algorithm,
  max_time,
  n_visits,
  amp_noise,
  period_noise,
  lpme_error,
  pme_error,
  principal_curve_error,
  pme_ssd,
  principal_curve_ssd,
  data_error,
  time_trend
)
error_df2 <- error_df %>%
  dplyr::select(-pme_ssd, -principal_curve_ssd) %>%
  pivot_longer(
  lpme_error:data_error,
  names_to = "model_type",
  values_to = "error"
) %>%
mutate(
  amp_noise_fct = as.factor(amp_noise),
  period_noise_fct = as.factor(period_noise),
  model_type = ifelse(
    model_type == "lpme_error",
    "LPME",
    ifelse(
      model_type == "pme_error",
      "PME",
      ifelse(
        model_type == "principal_curve_error",
        "PC",
        ifelse(
          model_type == "PME_error",
          "PME2",
          "Data"
        )
      )
    )
  )
) %>%
rename(Model = model_type) %>%
  mutate(
    amp_noise = as.factor(amp_noise),
    period_noise = as.factor(period_noise)
  )

error_df <- error_df %>%
  mutate(diff_error = pme_error - lpme_error)

p <- error_df2 %>%
  filter(sim_case == 7) %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(amp_noise), y = log(error), fill = Model)) +
  xlab("Amplitude Noise") +
  ylab("Mean Squared Distance (Log)") +
  ggtitle("LPME vs. PME, Simulation Mean Squared Distance")

ggsave("simulations_amp_error.png", plot = p)

p2 <- ggplot(error_df2) +
  geom_boxplot(aes(x = as.factor(period_noise), y = log(error), fill = Model)) +
  xlab("Shape Noise") +
  ylab("Mean Squared Distance (Log)")

# ggsave("partial_sim_results.png", p)

# error_list[[1]]$plot
# error_list[[6]]$plot

error_df2 %>%
  group_by(sim_case, time_trend, Model) %>%
  summarize(
    mean_error = mean(error),
    sd_error = sd(error)
  )

error_df2 %>%
  group_by(sim_case, time_trend, Model) %>%
  summarize(
    median_error = median(error),
    iqr_error = IQR(error)
  )

