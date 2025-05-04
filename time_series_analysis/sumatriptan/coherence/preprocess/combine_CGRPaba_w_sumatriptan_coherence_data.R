library(dplyr)


cgrpaba_coherence_df_filename = paste0(
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

sumatriptan_coherence_df_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/sumatriptan",
  "/coherence",
  "/data",
  "/sumatriptan_normalized_mean_logodds_coherence_combined_bla_cea_to_amy.RData");

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
  "/CGRPaba_w_sumatriptan_normalized_mean_logodds_coherence_bla_cea_to_amy.RData");

cgrpaba_coherence_ts_df = tibble(readRDS(cgrpaba_coherence_df_filename) );
sumatriptan_coherence_ts_df = tibble(readRDS(sumatriptan_coherence_df_filename) );

stopifnot(
  ncol(cgrpaba_coherence_ts_df$logodds_coherence_baseline) ==
  ncol(sumatriptan_coherence_ts_df$logodds_coherence_baseline) );

stopifnot(
  ncol(cgrpaba_coherence_ts_df$logodds_coherence_post_inj) ==
  ncol(sumatriptan_coherence_ts_df$logodds_coherence_post_inj) );

stopifnot(length(setdiff(names(cgrpaba_coherence_ts_df), names(sumatriptan_coherence_ts_df) )) == 0);
stopifnot(length(setdiff(names(sumatriptan_coherence_ts_df), names(cgrpaba_coherence_ts_df) )) == 0);

n_cgrp_days = length(unique(cgrpaba_coherence_ts_df$day) );
sumatriptan_coherence_ts_df$day = sumatriptan_coherence_ts_df$day + n_cgrp_days;

cgrp_mouse_ids = unique(cgrpaba_coherence_ts_df$mouse_id);
sumatriptan_mouse_ids = unique(sumatriptan_coherence_ts_df$mouse_id);

if (length(setdiff(sumatriptan_mouse_ids, cgrp_mouse_ids) ) > 0) {
  warning("Some mouse IDs in the Sumatriptan dataset do not appear in the CGRPaba dataset.");
}

cgrpaba_coherence_ts_df = filter(
  cgrpaba_coherence_ts_df,
  is.element(mouse_id, sumatriptan_mouse_ids) );

sumatriptan_coherence_ts_df = bind_rows(cgrpaba_coherence_ts_df, sumatriptan_coherence_ts_df);

saveRDS(sumatriptan_coherence_ts_df, save_filename);
