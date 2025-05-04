library(dplyr)
library(ggplot2)
library(purrr)
library(reshape2)
library(stringr)
library(tidyfun)
library(tidyr)
library(zoo)


coherence_ts_df_filename = paste0(
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

plots_save_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/CGRPaba",
  "/coherence",
  "/rolling_means_analysis",
  "/rolling_means_plots/");

if (!dir.exists(plots_save_dir) ) {
  dir.create(plots_save_dir, recursive=TRUE);

  stopifnot(dir.exists(plots_save_dir) );
}

coherence_ts_df = tibble(readRDS(coherence_ts_df_filename) );

mice_to_leave_out = c("MouseM36", "MouseMYF03");

coherence_ts_df = filter(coherence_ts_df, !is.element(mouse_id, mice_to_leave_out) );

time_baseline = ncol(coherence_ts_df$logodds_coherence_baseline);
time_post_inj = ncol(coherence_ts_df$logodds_coherence_post_inj);
post_inj_ts_args = ((1:time_post_inj) + time_baseline) / 60;

na_col_ixs = as.numeric(
  which(
    apply(
      coherence_ts_df$logodds_coherence_post_inj,
      MARGIN=2,
      FUN=\(col_ix) any(is.na(col_ix) ))));

if (length(na_col_ixs) > 0) {
  post_inj_ts_args = post_inj_ts_args[-na_col_ixs];
  coherence_ts_df$logodds_coherence_post_inj = coherence_ts_df$logodds_coherence_post_inj[,-na_col_ixs];
}

coherence_ts_df$logodds_coherence_post_inj = tfd(
  coherence_ts_df$logodds_coherence_post_inj,
  post_inj_ts_args);

ylim = c(-0.5, 0.8);
x_axis_ticks = seq(
  floor(min(post_inj_ts_args) ),
  floor(max(post_inj_ts_args) ),
  by=10);

reg_freq_group_dfs = (
  coherence_ts_df %>%
  group_by(region_pair, freq_band) %>%
  nest() );

n_reg_freqs = nrow(reg_freq_group_dfs);
window_size = 30;

for (rf_ix in 1:n_reg_freqs) {
  reg_ix = reg_freq_group_dfs$region_pair[rf_ix];
  freq_band_ix = reg_freq_group_dfs$freq_band[rf_ix];

  grp_df_ix = (
    reg_freq_group_dfs$data[[rf_ix]] %>%
    select(mouse_id, exp, logodds_coherence_post_inj) %>%
    mutate(
      smoothed_post_inj=tf_smooth(
        logodds_coherence_post_inj,
        method="rollmean",
        k=window_size,
        align="center") ));

  grp_df_ix$exp = as.factor(grp_df_ix$exp);

  cgrp_comp_df_ix = (
    bind_rows(
      grp_df_ix %>%
      select(mouse_id, exp, smoothed_post_inj) %>%
      filter(exp == "cgrp"),
      grp_df_ix %>%
      filter(exp == "veh") %>%
      group_by(mouse_id, exp) %>%
      summarize(
        smoothed_post_inj=mean(smoothed_post_inj, na.rm=TRUE),
        .groups="drop") ) %>%
    group_by(exp) %>%
    summarize(
      mean_post_inj=mean(smoothed_post_inj, na.rm=TRUE),
      se_post_inj=sd(smoothed_post_inj, na.rm=TRUE) / sqrt(n() ),
      .groups="drop") %>%
    mutate(
      lb_mean=mean_post_inj - se_post_inj,
      ub_mean=mean_post_inj + se_post_inj) %>%
    rename(treatment=exp) );

  cgrp_comp_plot_ix =  (
    cgrp_comp_df_ix %>%
    ggplot(aes(y=mean_post_inj, color=treatment) ) +
    coord_cartesian(ylim=ylim) +
    scale_x_continuous(breaks=x_axis_ticks) +
    theme(
      axis.text.x=element_text(size=28),
      axis.text.y=element_text(size=28),
      axis.title=element_text(size=20),
      legend.text=element_text(size=20),
      legend.title=element_text(size=20),
      plot.title=element_text(size=20) ) +
    geom_spaghetti() +
    geom_errorband(aes(ymax=ub_mean, ymin=lb_mean, fill=treatment) ) +
    labs(
      x="Time (minutes)",
      y="Mean log-odds coherence",
      title=paste(
        reg_ix,
        freq_band_ix,
        window_size,
        "second rolling mean") ));

  save_filename = paste0(
    plots_save_dir,
    sprintf("/%s_%d.png", reg_ix, freq_band_ix) );

  ggsave(
    save_filename,
    cgrp_comp_plot_ix,
    width=14,
    height=9);
}
