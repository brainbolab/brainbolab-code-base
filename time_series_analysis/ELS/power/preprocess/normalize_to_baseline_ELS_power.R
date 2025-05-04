r_data_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/ELS",
  "/power",
  "/data",
  "/ELS_mean_logpower_combined_bla_cea_to_amy.RData");

save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/ELS",
  "/power",
  "/data",
  "/ELS_normalized_mean_logpower_combined_bla_cea_to_amy.RData");

r_data = readRDS(r_data_filename);

mean_baseline = rowMeans(r_data$logpower_baseline, na.rm=TRUE);
r_data$logpower_baseline = (
  r_data$logpower_baseline -
  matrix(
    rep(mean_baseline, ncol(r_data$logpower_baseline) ),
    nrow=nrow(r_data$logpower_baseline) ));

r_data$logpower_post_inj = (
  r_data$logpower_post_inj -
  matrix(
    rep(mean_baseline, ncol(r_data$logpower_post_inj) ),
    nrow=nrow(r_data$logpower_post_inj) ));

saveRDS(r_data, save_filename);
