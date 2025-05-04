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
  "/ELS",
  "/coherence",
  "/data",
  "/ELS_normalized_mean_logodds_coherence_combined_bla_cea_to_amy.RData");

plots_save_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/ELS",
  "/coherence",
  "/gam_analysis",
  "/gam_plots/");

if (!dir.exists(plots_save_dir) ) {
  dir.create(plots_save_dir, recursive=TRUE);

  stopifnot(dir.exists(plots_save_dir) );
}

els_vs_control_effect_plots_save_dir = paste0(
  plots_save_dir,
  "/ELS_effect_plots/");

if (!dir.exists(els_vs_control_effect_plots_save_dir) ) {
  dir.create(els_vs_control_effect_plots_save_dir);

  stopifnot(dir.exists(els_vs_control_effect_plots_save_dir) );
}

els_cgrp_interaction_effect_plots_save_dir = paste0(
  plots_save_dir,
  "/ELS_CGRP_interaction_effect_plots/");

if (!dir.exists(els_cgrp_interaction_effect_plots_save_dir) ) {
  dir.create(els_cgrp_interaction_effect_plots_save_dir);

  stopifnot(dir.exists(els_cgrp_interaction_effect_plots_save_dir) );
}

effects_save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/ELS",
  "/coherence",
  "/gam_analysis",
  "/ELS_GAM_analysis_effects.csv");

pvalues_save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/ELS",
  "/coherence",
  "/gam_analysis",
  "/ELS_GAM_analysis_pvalues.csv");

all_regions_fitted_save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/ELS",
  "/coherence",
  "/gam_analysis",
  "/ELS_GAM_analysis_fitted_w_se.RData");

coherence_ts_df = tibble(readRDS(coherence_ts_df_filename) );

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

time_cutoff = 75; # In minutes.
t_dim = sum(plot_post_inj_ts_args <= time_cutoff);

plot_post_inj_ts_args = plot_post_inj_ts_args[1:t_dim];
post_inj_ts_args = post_inj_ts_args[1:t_dim];
coherence_ts_df$logodds_coherence_post_inj = coherence_ts_df$logodds_coherence_post_inj[,1:t_dim];

coherence_ts_df$logodds_coherence_post_inj = tfd(
  coherence_ts_df$logodds_coherence_post_inj,
  post_inj_ts_args);

coherence_ts_df$mouse_id = as.factor(coherence_ts_df$mouse_id);
coherence_ts_df$cgrp = factor(coherence_ts_df$exp, levels=c("veh", "cgrp"), ordered=TRUE);
coherence_ts_df$els_w_cgrp = factor(0, levels=c(0, 1), ordered=TRUE);
els_cgrp_intn_ixs = (coherence_ts_df$els == "els") & (coherence_ts_df$exp == "cgrp");
coherence_ts_df$els_w_cgrp[els_cgrp_intn_ixs] = 1;
coherence_ts_df$els = factor(coherence_ts_df$els, levels=c("control", "els"), ordered=TRUE);

# Group the data by region/frequency band.
reg_freq_group_dfs = (
  coherence_ts_df %>%
  group_by(region_pair, freq_band) %>%
  nest() %>%
  mutate(
    els_par_pval=NA,
    els_spval=NA,
    els_cgrp_intn_par_pval=NA,
    els_cgrp_intn_spval=NA,
    els_effect=NA,
    els_cgrp_intn_effect=NA) );

n_reg_freqs = nrow(reg_freq_group_dfs);

reg_freq_group_dfs$els_fitted_vals = matrix(NA, n_reg_freqs, t_dim);
reg_freq_group_dfs$els_se = matrix(NA, n_reg_freqs, t_dim);
reg_freq_group_dfs$els_cgrp_intn_fitted_vals = matrix(NA, n_reg_freqs, t_dim);
reg_freq_group_dfs$els_cgrp_intn_se = matrix(NA, n_reg_freqs, t_dim);

# These variables pertain to plot customization.
ylim = c(-0.8, 0.6);
x_axis_ticks = seq(
  floor(min(plot_post_inj_ts_args) ),
  floor(max(plot_post_inj_ts_args) ),
  by=10);

for (rf_ix in 1:n_reg_freqs) {
  reg_ix = reg_freq_group_dfs$region_pair[rf_ix];
  freq_band_ix = reg_freq_group_dfs$freq_band[rf_ix];

  grp_df_ix = reg_freq_group_dfs$data[[rf_ix]];

  els_vs_control_gam_df = (
    grp_df_ix %>%
    select(mouse_id, logodds_coherence_post_inj, els, cgrp, els_w_cgrp) %>%
    tf_unnest(cols="logodds_coherence_post_inj") %>%
    rename(time=logodds_coherence_post_inj_arg, reg_freq=logodds_coherence_post_inj_value) );

  els_vs_control_gam_formula = as.formula(
    paste0(
      "reg_freq ~ ",
      "s(mouse_id, bs='re') + ",
      "s(mouse_id, time, bs='re') + ",
      "s(time) + ",
      "els + ",
      "cgrp + ",
      "els_w_cgrp + ",
      "s(time, by=els) + ",
      "s(time, by=cgrp) + ",
      "s(time, by=els_w_cgrp)") );

  els_vs_control_gam_mod = gam(els_vs_control_gam_formula, data=els_vs_control_gam_df);

  reg_freq_group_dfs$els_par_pval[rf_ix] = max(
    summary(els_vs_control_gam_mod)$p.pv[["els.L"]],
    2e-16);

  reg_freq_group_dfs$els_spval[rf_ix] = max(
    summary(els_vs_control_gam_mod)$s.pv[4],
    2e-16);

  reg_freq_group_dfs$els_cgrp_intn_par_pval[rf_ix] = max(
    summary(els_vs_control_gam_mod)$p.pv[["els_w_cgrp.L"]],
    2e-16);

  reg_freq_group_dfs$els_cgrp_intn_spval[rf_ix] = max(
    summary(els_vs_control_gam_mod)$s.pv[6],
    2e-16);

  reg_freq_group_dfs$els_effect[rf_ix] = els_vs_control_gam_mod$coefficients[["els.L"]];
  reg_freq_group_dfs$els_cgrp_intn_effect[rf_ix] = els_vs_control_gam_mod$coefficients[["els_w_cgrp.L"]];

  cat(
    sprintf(
      paste0(
        "P-values for %s %d:\n",
        "          ELS par pv = %f\n",
        "       ELS smooth pv = %f\n",
        "     ELS:CGRP par pv = %f\n",
        "  ELS:CGRP smooth pv = %f\n"),
      reg_ix,
      freq_band_ix,
      reg_freq_group_dfs$els_par_pval[rf_ix],
      reg_freq_group_dfs$els_spval[rf_ix],
      reg_freq_group_dfs$els_cgrp_intn_par_pval[rf_ix],
      reg_freq_group_dfs$els_cgrp_intn_spval[rf_ix]) );

  els_vs_control_gam_fit = predict(
    els_vs_control_gam_mod,
    data.frame(
      time=post_inj_ts_args,
      mouse_id=grp_df_ix$mouse_id[1],
      els="els",
      cgrp="cgrp",
      els_w_cgrp=1),
    type="terms",
    se.fit=TRUE);

  reg_freq_group_dfs$els_fitted_vals[rf_ix,] = els_vs_control_gam_fit$fit[,"s(time):elsels"];
  reg_freq_group_dfs$els_se[rf_ix,] = els_vs_control_gam_fit$se.fit[,"s(time):elsels"];
  reg_freq_group_dfs$els_cgrp_intn_fitted_vals[rf_ix,] = els_vs_control_gam_fit$fit[,"s(time):els_w_cgrp1"];
  reg_freq_group_dfs$els_cgrp_intn_se[rf_ix,] = els_vs_control_gam_fit$se.fit[,"s(time):els_w_cgrp1"];

  gam_fit_plot_df = tibble(NA);
  gam_fit_plot_df$gam_fit = tfd(
    reg_freq_group_dfs$els_fitted_vals[rf_ix,],
    plot_post_inj_ts_args);

  gam_fit_plot_df$gam_fit_ub = tfd(
    reg_freq_group_dfs$els_fitted_vals[rf_ix,] + (2 * reg_freq_group_dfs$els_se[rf_ix,]),
    plot_post_inj_ts_args);

  gam_fit_plot_df$gam_fit_lb = tfd(
    reg_freq_group_dfs$els_fitted_vals[rf_ix,] - (2 * reg_freq_group_dfs$els_se[rf_ix,]),
    plot_post_inj_ts_args);

  els_vs_control_plot_ix = (
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
      y="ELS - control effect mean log-odds coherence",
      title=sprintf(
        "GAM Analysis of ELS Effect on %s %d",
        reg_ix,
        freq_band_ix) ));

  save_filename = paste0(
    els_vs_control_effect_plots_save_dir,
    sprintf("/%s_%d.png", reg_ix, freq_band_ix) );

  ggsave(
    save_filename,
    els_vs_control_plot_ix,
    width=13,
    height=9);

  gam_fit_plot_df = tibble(NA);
  gam_fit_plot_df$gam_fit = tfd(
    reg_freq_group_dfs$els_cgrp_intn_fitted_vals[rf_ix,],
    plot_post_inj_ts_args);

  gam_fit_plot_df$gam_fit_ub = tfd(
    reg_freq_group_dfs$els_cgrp_intn_fitted_vals[rf_ix,] + (2 * reg_freq_group_dfs$els_cgrp_intn_se[rf_ix,]),
    plot_post_inj_ts_args);

  gam_fit_plot_df$gam_fit_lb = tfd(
    reg_freq_group_dfs$els_cgrp_intn_fitted_vals[rf_ix,] - (2 * reg_freq_group_dfs$els_cgrp_intn_se[rf_ix,]),
    plot_post_inj_ts_args);

  els_cgrp_intn_plot_ix = (
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
      y="ELS:CGRP interaction effect mean log-odds coherence",
      title=sprintf(
        "GAM Analysis of ELS:CGRP Interaction Effect on %s %d",
        reg_ix,
        freq_band_ix) ));

  save_filename = paste0(
    els_cgrp_interaction_effect_plots_save_dir,
    sprintf("/%s_%d.png", reg_ix, freq_band_ix) );

  ggsave(
    save_filename,
    els_cgrp_intn_plot_ix,
    width=13,
    height=9);
}

write.csv(
  reg_freq_group_dfs %>%
  select(region_pair, freq_band, els_effect, els_cgrp_intn_effect),
  effects_save_filename,
  row.names=FALSE);

write.csv(
  reg_freq_group_dfs %>%
  select(
    region_pair,
    freq_band,
    els_par_pval,
    els_spval,
    els_cgrp_intn_par_pval,
    els_cgrp_intn_spval),
  pvalues_save_filename,
  row.names=FALSE);

saveRDS(
  reg_freq_group_dfs %>%
  select(
    region_pair,
    freq_band,
    els_fitted_vals,
    els_se,
    els_cgrp_intn_fitted_vals,
    els_cgrp_intn_se),
  all_regions_fitted_save_filename);
