# Run this script after running 'concatenate_separate_power_ts.m' which produces the CSV files needed for
# this script.

library(dplyr)


# Directory where the CSV files constructed from 'concatenate_separate_power_ts.m' are located.
data_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/average_power",
  "/ELS_all_data",
  "/CSV_data/");

# Save filename for the R dataframe created by combining the CSV files.
save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/average_power",
  "/ELS_all_data",
  "/R_data",
  "/ELS_mean_logpower_combined.RData");

# Names of the CSV files to combine into an R dataframe.
data_filenames = c(
  "Day1_ELS_power_fused_lasso_analysis.csv",
  "Day2_ELS_power_fused_lasso_analysis.csv",
  "Day3_ELS_power_fused_lasso_analysis.csv",
  "Day1_ELS_control_power_fused_lasso_analysis.csv",
  "Day2_ELS_control_power_fused_lasso_analysis.csv",
  "Day3_ELS_control_power_fused_lasso_analysis.csv");

# Labels that should be added to each row of the resulting combined data. These must match the order of
# the CSV data filenames.
exp_labels = c("veh", "cgrp", "veh", "veh", "cgrp", "veh");
day_labels = c(1:3, 1:3);
els_labels = c("els", "els", "els", "control", "control", "control");

# It's a good idea to look at the 'filenames_info_table' object to make sure that the labels match to
# the filenames.
filenames_info_table = data.frame(
  filename=data_filenames,
  exp=exp_labels,
  day=day_labels,
  els=els_labels);

n_filenames = nrow(filenames_info_table);

# The variables in the 'filenames_info_table' dataframe that should be added to each row of the combined data.
addl_tags = c("exp", "day", "els");
n_tags = length(addl_tags);

# Each phase of the time series' that should be separated into its own matrix, along with the start and end
# indices.
exp_phases_info = data.frame(
  exp_phase=c("baseline", "post_inj"),
  ix_start=c(1, 1802),
  ix_end=c(1801, 7204) );

n_exp_phases = nrow(exp_phases_info);

expected_nonpow_colnames = c(
  "mouse_id",
  "freq_band",
  "region");

rdata_combined = NULL;

for (fx in 1:n_filenames) {
  filename_fx = paste0(data_dir, "/", filenames_info_table$filename[fx]);
  data_fx = read.csv(filename_fx, header=TRUE);
  nonpow_col_ixs = which(is.element(colnames(data_fx), expected_nonpow_colnames) );

  stopifnot(length(nonpow_col_ixs) == length(expected_nonpow_colnames) );

  rdata_fx = data_fx[,nonpow_col_ixs];

  for (tx in 1:n_tags) {
    tag_tx = addl_tags[tx];
    rdata_fx[[tag_tx]] = filenames_info_table[[tag_tx]][fx];
  }

  logpower_mat = as.matrix(data_fx[,-nonpow_col_ixs]);
  colnames(logpower_mat) = NULL;

  for (px in 1:n_exp_phases) {
    varname_px = paste0("logpower_", exp_phases_info$exp_phase[px]);
    exp_phase_px_ixs = exp_phases_info$ix_start[px]:exp_phases_info$ix_end[px];
    rdata_fx[[varname_px]] = logpower_mat[,exp_phase_px_ixs];
  }

  if (is.null(rdata_combined) ) {
    rdata_combined = rdata_fx;
  } else {
    rdata_combined = bind_rows(rdata_combined, rdata_fx);
  }
}

saveRDS(rdata_combined, save_filename);
