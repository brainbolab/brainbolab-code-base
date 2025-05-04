library(tidyverse)


r_data_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/sumatriptan",
  "/coherence",
  "/data",
  "/sumatriptan_mean_coherence_combined.RData");

save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/sumatriptan",
  "/coherence",
  "/data",
  "/sumatriptan_mean_coherence_combined_bla_cea_to_amy.RData");

r_data = readRDS(r_data_filename);

combine_regs_pattern = "(BLA)|(CeA)";
link_regs_pattern = "(Acc)|(LPBN)|(MD_thal)|(Po)|(RPBN)|(VPM)";

amy_groups = (
  r_data %>%
  filter(str_detect(region_pair, combine_regs_pattern) ) %>%
  mutate(link_region=str_extract(region_pair, link_regs_pattern) ) %>%
  filter(!is.na(link_region) ) %>%
  group_by(
    mouse_id,
    freq_band,
    exp,
    day,
    link_region) %>%
  nest() %>%
  mutate(region_pair=NA) %>%
  as.data.frame() );

n_obs = nrow(amy_groups);

amy_groups$coherence_baseline = matrix(NA, n_obs, dim(r_data$coherence_baseline)[2]);
amy_groups$coherence_post_inj = matrix(NA, n_obs, dim(r_data$coherence_post_inj)[2]);

for (ix in 1:n_obs) {
  amy_groups$coherence_baseline[ix,] = colMeans(amy_groups$data[[ix]]$coherence_baseline);
  amy_groups$coherence_post_inj[ix,] = colMeans(amy_groups$data[[ix]]$coherence_post_inj);

  amy_groups$region_pair[ix] = paste0(
    sort(c("AMY", amy_groups$link_region[ix]) ),
    collapse="_x_");
}

r_data = bind_rows(
  filter(r_data, !str_detect(region_pair, combine_regs_pattern) ),
  select(amy_groups, !c(data, link_region) ));

saveRDS(r_data, save_filename);
