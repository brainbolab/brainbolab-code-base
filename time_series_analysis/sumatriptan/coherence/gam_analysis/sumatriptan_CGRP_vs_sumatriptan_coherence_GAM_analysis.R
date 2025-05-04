library(dplyr)
library(ggplot2)
library(mgcv)
library(reshape2)
library(tidyfun)
library(tidyr)


coherence_ts_df_filename = paste0(
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

plots_save_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/sumatriptan",
  "/coherence",
  "/gam_analysis",
  "/gam_plots/");

if (!dir.exists(plots_save_dir) ) {
  dir.create(plots_save_dir);

  stopifnot(dir.exists(plots_save_dir) );
}

pvalues_save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/sumatriptan",
  "/coherence",
  "/gam_analysis",
  "/sumatriptan_GAM_analysis_pvalues.csv");

all_regions_fitted_save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/sumatriptan",
  "/coherence",
  "/gam_analysis",
  "/sumatriptan_GAM_analysis_fitted_w_se.RData");

coherence_ts_df = tibble(readRDS(coherence_ts_df_filename) );

mice_to_leave_out = c("MouseM36", "MouseMYF03");

coherence_ts_df = filter(coherence_ts_df, !is.element(mouse_id, mice_to_leave_out) );

time_baseline = ncol(coherence_ts_df$logodds_coherence_baseline);
time_post_inj = ncol(coherence_ts_df$logodds_coherence_post_inj);
plot_post_inj_ts_args = ((1:time_post_inj) + time_baseline) / 60;
post_inj_ts_args = scale(0:(time_post_inj - 1) );

na_col_ixs = as.numeric(
  which(
    apply(
      coherence_ts_df$logodds_coherence_post_inj,
      MARGIN=2,
      FUN=\(col_ix) any(is.na(col_ix) ))));

if (length(na_col_ixs) > 0) {
  plot_post_inj_ts_args = plot_post_inj_ts_args[-na_col_ixs];
  post_inj_ts_args = post_inj_ts_args[-na_col_ixs];
  coherence_ts_df$logodds_coherence_post_inj = coherence_ts_df$logodds_coherence_post_inj[,-na_col_ixs];
}

time_cutoff = 50; # In minutes.
t_dim = sum(plot_post_inj_ts_args <= time_cutoff);

plot_post_inj_ts_args = plot_post_inj_ts_args[1:t_dim];
post_inj_ts_args = post_inj_ts_args[1:t_dim];
coherence_ts_df$logodds_coherence_post_inj = coherence_ts_df$logodds_coherence_post_inj[,1:t_dim];

coherence_ts_df$logodds_coherence_post_inj = tfd(
  coherence_ts_df$logodds_coherence_post_inj,
  post_inj_ts_args);

coherence_ts_df$mouse_id = as.factor(coherence_ts_df$mouse_id);

# Group the data by region_pair/frequency band.
reg_freq_group_dfs = (
  coherence_ts_df %>%
  group_by(region_pair, freq_band) %>%
  nest() %>%
  mutate(s_pval=NA, int_pval=NA) );

n_reg_freqs = nrow(reg_freq_group_dfs);

reg_freq_group_dfs$fitted_vals = matrix(NA, n_reg_freqs, t_dim);
reg_freq_group_dfs$se = matrix(NA, n_reg_freqs, t_dim);

# These variables pertain to plot customization.
ylim = c(-0.5, 0.4);
x_axis_ticks = seq(
  floor(min(plot_post_inj_ts_args) ),
  floor(max(plot_post_inj_ts_args) ),
  by=10);

for (rf_ix in 1:n_reg_freqs) {
  reg_ix = reg_freq_group_dfs$region_pair[rf_ix];
  freq_band_ix = reg_freq_group_dfs$freq_band[rf_ix];

  grp_df_ix = reg_freq_group_dfs$data[[rf_ix]];

  # For each mouse, subtract the sumatriptan w/ CGRP treatment post injection time series
  # from its CGRP only treatment post injection time series.
  grp_diff_df_ix = (
    inner_join(
      grp_df_ix %>%
      filter(exp == "sumatriptan_w_cgrp") %>%
      select(mouse_id, logodds_coherence_post_inj) %>%
      rename(sumatriptan_w_cgrp_post_inj=logodds_coherence_post_inj),
      grp_df_ix %>%
      filter(exp == "cgrp") %>%
      select(mouse_id, logodds_coherence_post_inj) %>%
      rename(cgrp_post_inj=logodds_coherence_post_inj),
      by="mouse_id") %>%
    mutate(cgrp_minus_sumatriptan=cgrp_post_inj - sumatriptan_w_cgrp_post_inj) %>%
    select(mouse_id, cgrp_minus_sumatriptan) );

  # Melt resulting dataframe.
  cgrp_minus_sumatriptan_gam_df = (
    grp_diff_df_ix %>%
    tf_unnest(cols="cgrp_minus_sumatriptan") %>%
    rename(
      time=cgrp_minus_sumatriptan_arg,
      reg_freq_diff=cgrp_minus_sumatriptan_value) );

  cgrp_minus_sumatriptan_gam_mod = gam(
    reg_freq_diff ~ s(time),
    data=cgrp_minus_sumatriptan_gam_df);

  reg_freq_group_dfs$s_pval[rf_ix] = max(
    summary(cgrp_minus_sumatriptan_gam_mod)$s.pv,
    2e-16);

  reg_freq_group_dfs$int_pval[rf_ix] = max(
    summary(cgrp_minus_sumatriptan_gam_mod)$p.pv,
    2e-16);

  print(
    sprintf(
      "P-values for %s %d from mod w/o rand eff: s_pv = %f, int_pv = %f",
      reg_ix,
      freq_band_ix,
      summary(cgrp_minus_sumatriptan_gam_mod)$s.pv,
      summary(cgrp_minus_sumatriptan_gam_mod)$p.pv) );

  c_minus_s_gam_fit = predict(
    cgrp_minus_sumatriptan_gam_mod,
    data.frame(time=post_inj_ts_args),
    se.fit=TRUE);

  reg_freq_group_dfs$fitted_vals[rf_ix,] = c_minus_s_gam_fit$fit;
  reg_freq_group_dfs$se[rf_ix,] = c_minus_s_gam_fit$se.fit;

  gam_fit_plot_df = tibble(NA);
  gam_fit_plot_df$gam_fit = tfd(c_minus_s_gam_fit$fit, arg=plot_post_inj_ts_args);
  gam_fit_plot_df$gam_fit_ub = tfd(
    c_minus_s_gam_fit$fit + (2 * c_minus_s_gam_fit$se.fit),
    arg=plot_post_inj_ts_args);

  gam_fit_plot_df$gam_fit_lb = tfd(
    c_minus_s_gam_fit$fit - (2 * c_minus_s_gam_fit$se.fit),
    arg=plot_post_inj_ts_args);

  sumatriptan_vs_cgrp_effect_plot = (
    gam_fit_plot_df %>%
    ggplot(aes(y=gam_fit) ) +
    coord_cartesian(ylim=ylim) +
    scale_x_continuous(breaks=x_axis_ticks) +
    theme(
      axis.text.x=element_text(size=28),
      axis.text.y=element_text(size=28),
      axis.title=element_text(size=20),
      plot.title=element_text(size=20) ) +
    geom_spaghetti() +
    geom_errorband(aes(ymax=gam_fit_ub, ymin=gam_fit_lb) ) +
    geom_hline(
      yintercept=0,
      linetype="dashed",
      color="red") +
    labs(
      x="Time (minutes)",
      y="CGRP only - Sumatriptan w/ CGRP mean log-odds coherence",
      title=sprintf(
        "GAM Analysis of CGRP only vs Sumatriptan w/ CGRP effect on %s %d",
        reg_ix,
        freq_band_ix) ));

  save_filename = paste0(
    plots_save_dir,
    sprintf("/%s_%d.png", reg_ix, freq_band_ix) );

  ggsave(
    save_filename,
    sumatriptan_vs_cgrp_effect_plot,
    width=13,
    height=9);
}

write.csv(
  select(reg_freq_group_dfs, region_pair, freq_band, s_pval, int_pval),
  pvalues_save_filename,
  row.names=FALSE);

saveRDS(
  reg_freq_group_dfs %>%
  select(region_pair, freq_band, fitted_vals, se),
  all_regions_fitted_save_filename);
