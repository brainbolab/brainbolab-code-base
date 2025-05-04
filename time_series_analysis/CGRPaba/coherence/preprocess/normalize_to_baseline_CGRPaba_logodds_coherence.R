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
  "/CGRPaba_mean_logodds_coherence_combined_bla_cea_to_amy.RData");

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
  "/CGRPaba_normalized_mean_logodds_coherence_combined_bla_cea_to_amy.RData");

r_data = readRDS(r_data_filename);

mean_baseline = rowMeans(r_data$logodds_coherence_baseline, na.rm=TRUE);
r_data$logodds_coherence_baseline = (
  r_data$logodds_coherence_baseline -
  matrix(
    rep(mean_baseline, ncol(r_data$logodds_coherence_baseline) ),
    nrow=nrow(r_data$logodds_coherence_baseline) ));

r_data$logodds_coherence_post_inj = (
  r_data$logodds_coherence_post_inj -
  matrix(
    rep(mean_baseline, ncol(r_data$logodds_coherence_post_inj) ),
    nrow=nrow(r_data$logodds_coherence_post_inj) ));

saveRDS(r_data, save_filename);
