library(dplyr)


cgrpaba_power_df_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/CGRPaba",
  "/power",
  "/data",
  "/CGRPaba_nomalized_mean_logpower_combined_bla_cea_to_amy.RData");

sumatriptan_power_df_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/sumatriptan",
  "/power",
  "/data",
  "/sumatriptan_normalized_mean_logpower_combined_bla_cea_to_amy.RData");

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
  "/CGRPaba_w_sumatriptan_normalized_mean_logpower_bla_cea_to_amy.RData");

cgrpaba_power_ts_df = tibble(readRDS(cgrpaba_power_df_filename) );
sumatriptan_power_ts_df = tibble(readRDS(sumatriptan_power_df_filename) );

stopifnot(ncol(cgrpaba_power_ts_df$logpower_baseline) == ncol(sumatriptan_power_ts_df$logpower_baseline) );
stopifnot(ncol(cgrpaba_power_ts_df$logpower_post_inj) == ncol(sumatriptan_power_ts_df$logpower_post_inj) );
stopifnot(length(setdiff(names(cgrpaba_power_ts_df), names(sumatriptan_power_ts_df) )) == 0);
stopifnot(length(setdiff(names(sumatriptan_power_ts_df), names(cgrpaba_power_ts_df) )) == 0);

n_cgrp_days = length(unique(cgrpaba_power_ts_df$day) );
sumatriptan_power_ts_df$day = sumatriptan_power_ts_df$day + n_cgrp_days;

cgrp_mouse_ids = unique(cgrpaba_power_ts_df$mouse_id);
sumatriptan_mouse_ids = unique(sumatriptan_power_ts_df$mouse_id);

if (length(setdiff(sumatriptan_mouse_ids, cgrp_mouse_ids) ) > 0) {
  warning("Some mouse IDs in the Sumatriptan dataset do not appear in the CGRPaba dataset.");
}

cgrpaba_power_ts_df = filter(
  cgrpaba_power_ts_df,
  is.element(mouse_id, sumatriptan_mouse_ids) );

sumatriptan_power_ts_df = bind_rows(cgrpaba_power_ts_df, sumatriptan_power_ts_df);

saveRDS(sumatriptan_power_ts_df, save_filename);
