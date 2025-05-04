data_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/average_power",
  "/CGRPaba_data",
  "/fused_lasso_data240323/");

day_data_filenames = c(
  "mean_logpower_time_series_results_postinjcut_Day1_fused_lasso_analysis.csv",
  "mean_logpower_time_series_results_postinjcut_Day2_fused_lasso_analysis.csv",
  "mean_logpower_time_series_results_postinjcut_Day3_fused_lasso_analysis.csv");

baseline_range = c(1, 1801);

expected_nonpow_colnames = c("mouse_id", "freq_band", "region");

save_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/average_power",
  "/ELS_control_data",
  "/normalized_to_baseline_fused_lasso_data240516/");

if (!dir.exists(save_dir) ) {
  dir.create(save_dir);

  stopifnot(dir.exists(save_dir) );
}

n_days = length(day_data_filenames);

for (dx in 1:n_days) {
  next_day_data_fname = paste0(data_dir, day_data_filenames[dx]);
  next_day_data = read.csv(next_day_data_fname, header=TRUE);

  n_obs = nrow(next_day_data);

  nonpow_col_ixs = which(
    is.element(colnames(next_day_data), expected_nonpow_colnames) );

  stopifnot(length(nonpow_col_ixs) == length(expected_nonpow_colnames) );

  all_power_col_ixs = sort(
    setdiff(
      1:ncol(next_day_data),
      nonpow_col_ixs) );

  n_power = length(all_power_col_ixs);

  baseline_power_col_ixs = all_power_col_ixs[baseline_range[1]:baseline_range[2]];
  baseline_means = rowMeans(
    as.matrix(next_day_data[,baseline_power_col_ixs]),
    na.rm=TRUE);

  next_day_data[,all_power_col_ixs] = (
    next_day_data[,all_power_col_ixs] -
    matrix(rep(baseline_means, n_power), n_obs, n_power) );

  save_filename = paste0(save_dir, "normalized_to_baseline_", day_data_filenames[dx]);
  write.csv(next_day_data, save_filename, row.names=FALSE);
}
