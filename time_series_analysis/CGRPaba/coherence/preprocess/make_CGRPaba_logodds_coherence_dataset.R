r_data_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/CGRPaba",
  "/coherence",
  "/data",
  "/CGRPaba_mean_coherence_combined_bla_cea_to_amy.RData");

save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/CGRPaba",
  "/coherence",
  "/data",
  "/CGRPaba_mean_logodds_coherence_combined_bla_cea_to_amy.RData");

r_data = readRDS(r_data_filename);

r_data$logodds_coherence_baseline = (
  log(r_data$coherence_baseline) -
  log(1 - r_data$coherence_baseline) );

r_data$logodds_coherence_post_inj = (
  log(r_data$coherence_post_inj) -
  log(1 - r_data$coherence_post_inj) );

saveRDS(
  dplyr::select(r_data, !c(coherence_baseline, coherence_post_inj) ),
  save_filename);
