library(tidyverse)


r_data_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/sumatriptan",
  "/power",
  "/data",
  "/sumatriptan_mean_logpower_combined.RData");

save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/sumatriptan",
  "/power",
  "/data",
  "/sumatriptan_mean_logpower_combined_bla_cea_to_amy.RData");

r_data = readRDS(r_data_filename);

amy_groups = (
  r_data %>%
  filter(region == "BLA" | region == "CeA") %>%
  group_by(
    mouse_id,
    freq_band,
    exp,
    day) %>%
  nest() %>%
  mutate(region="AMY") );

n_obs = nrow(amy_groups);

amy_groups$logpower_baseline = matrix(NA, n_obs, dim(r_data$logpower_baseline)[2]);
amy_groups$logpower_post_inj = matrix(NA, n_obs, dim(r_data$logpower_post_inj)[2]);

for (ix in 1:n_obs) {
  amy_groups$logpower_baseline[ix,] = colMeans(amy_groups$data[[ix]]$logpower_baseline);
  amy_groups$logpower_post_inj[ix,] = colMeans(amy_groups$data[[ix]]$logpower_post_inj);
}

r_data = bind_rows(
  filter(r_data, region != "BLA" & region != "CeA"),
  select(amy_groups, !data) );

saveRDS(r_data, save_filename);
