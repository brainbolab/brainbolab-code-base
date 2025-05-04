library(dplyr)
library(ggplot2)
library(tidyfun)
library(tidyr)


power_ts_df_filename = paste0(
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

plots_save_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/ELS",
  "/power",
  "/rolling_means_analysis",
  "/rolling_means_plots/");

if (!dir.exists(plots_save_dir) ) {
  dir.create(plots_save_dir, recursive=TRUE);

  stopifnot(dir.exists(plots_save_dir) );
}

power_ts_df = tibble(readRDS(power_ts_df_filename) );

time_baseline = ncol(power_ts_df$logpower_baseline);
time_post_inj = ncol(power_ts_df$logpower_post_inj);
post_inj_ts_args = ((1:time_post_inj) + time_baseline) / 60;

na_col_ixs = as.numeric(
  which(
    apply(
      power_ts_df$logpower_post_inj,
      MARGIN=2,
      FUN=\(col_ix) any(is.na(col_ix) ))));

if (length(na_col_ixs) > 0) {
  post_inj_ts_args = post_inj_ts_args[-na_col_ixs];
  power_ts_df$logpower_post_inj = power_ts_df$logpower_post_inj[,-na_col_ixs];
}

time_cutoff = 75; # In minutes.
t_dim = sum(post_inj_ts_args <= time_cutoff);

post_inj_ts_args = post_inj_ts_args[1:t_dim];
power_ts_df$logpower_post_inj = power_ts_df$logpower_post_inj[,1:t_dim];

power_ts_df$logpower_post_inj = tfd(
  power_ts_df$logpower_post_inj,
  post_inj_ts_args);

ylim = c(-0.2, 0.15);
x_axis_ticks = seq(
  floor(min(post_inj_ts_args) ),
  floor(max(post_inj_ts_args) ),
  by=10);

reg_freq_group_dfs = (
  power_ts_df %>%
  group_by(region, freq_band) %>%
  nest() );

n_reg_freqs = nrow(reg_freq_group_dfs);
window_size = 30;

for (rf_ix in 1:n_reg_freqs) {
  reg_ix = reg_freq_group_dfs$region[rf_ix];
  freq_band_ix = reg_freq_group_dfs$freq_band[rf_ix];

  grp_df_ix = (
    reg_freq_group_dfs$data[[rf_ix]] %>%
    select(mouse_id, els, exp, logpower_post_inj) %>%
    mutate(
      smoothed_post_inj=tf_smooth(
        logpower_post_inj,
        method="rollmean",
        k=window_size,
        align="center") ));

  els_cgrp_comp_df_ix = (
    bind_rows(
      grp_df_ix %>%
      select(mouse_id, els, exp, smoothed_post_inj) %>%
      filter(exp == "cgrp"),
      grp_df_ix %>%
      filter(exp == "veh") %>%
      group_by(mouse_id, els, exp) %>%
      summarize(
        smoothed_post_inj=mean(smoothed_post_inj, na.rm=TRUE),
        .groups="drop") ) %>%
    group_by(els, exp) %>%
    summarize(
      mean_post_inj=mean(smoothed_post_inj, na.rm=TRUE),
      se_post_inj=sd(smoothed_post_inj, na.rm=TRUE) / sqrt(n() ),
      .groups="drop") %>%
    mutate(
      lb_mean=mean_post_inj - se_post_inj,
      ub_mean=mean_post_inj + se_post_inj) );

  els_cgrp_comp_df_ix$els[els_cgrp_comp_df_ix$els == "els"] = "ELS";
  els_cgrp_comp_df_ix$els[els_cgrp_comp_df_ix$els == "control"] = "Control";
  els_cgrp_comp_df_ix$exp[els_cgrp_comp_df_ix$exp == "cgrp"] = "CGRP";
  els_cgrp_comp_df_ix$exp[els_cgrp_comp_df_ix$exp == "veh"] = "Vehicle";
  els_cgrp_comp_df_ix = mutate(
    els_cgrp_comp_df_ix,
    treatment=as.factor(paste(els, exp) ));

  els_cgrp_comp_plot_ix =  (
    els_cgrp_comp_df_ix %>%
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
      y="Mean log-power",
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
    els_cgrp_comp_plot_ix,
    width=18,
    height=9);
}
