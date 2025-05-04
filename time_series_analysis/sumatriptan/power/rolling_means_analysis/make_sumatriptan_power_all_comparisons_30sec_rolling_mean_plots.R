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
  "/sumatriptan",
  "/power",
  "/data",
  "/CGRPaba_w_sumatriptan_normalized_mean_logpower_bla_cea_to_amy.RData");

plots_save_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/sumatriptan",
  "/power",
  "/rolling_means_analysis",
  "/rolling_means_plots/");

if (!dir.exists(plots_save_dir) ) {
  dir.create(plots_save_dir, recursive=TRUE);

  stopifnot(dir.exists(plots_save_dir) );
}

power_ts_df = tibble(readRDS(power_ts_df_filename) );

mice_to_leave_out = c("MouseM36", "MouseMYF03");

power_ts_df = filter(power_ts_df, !is.element(mouse_id, mice_to_leave_out) );

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

ylim = c(-0.15, 0.15);
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

line_colors = c("#F8766D", "#7CAE00", "#C77CFF", "#00BFC4");

for (rf_ix in 1:n_reg_freqs) {
  reg_ix = reg_freq_group_dfs$region[rf_ix];
  freq_band_ix = reg_freq_group_dfs$freq_band[rf_ix];

  grp_df_ix = (
    reg_freq_group_dfs$data[[rf_ix]] %>%
    select(mouse_id, exp, logpower_post_inj) %>%
    mutate(
      smoothed_post_inj=tf_smooth(
        logpower_post_inj,
        method="rollmean",
        k=window_size,
        align="center") ));

  grp_df_ix$exp = as.factor(grp_df_ix$exp);

  all_treatments_df_ix = (
    bind_rows(
      grp_df_ix %>%
      select(mouse_id, exp, smoothed_post_inj) %>%
      filter(exp != "veh"),
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

  all_treatments_df_ix$treatment = as.character(all_treatments_df_ix$treatment);
  all_treatments_df_ix$treatment[all_treatments_df_ix$treatment == "cgrp"] = "CGRP";
  all_treatments_df_ix$treatment[all_treatments_df_ix$treatment == "sumatriptan"] = "Sumatriptan";
  all_treatments_df_ix$treatment[all_treatments_df_ix$treatment == "sumatriptan_w_cgrp"] = "CGRP w Sumatriptan";
  all_treatments_df_ix$treatment[all_treatments_df_ix$treatment == "veh"] = "Vehicle";

  all_treatments_df_ix = all_treatments_df_ix[order(all_treatments_df_ix$treatment),];
  all_treatments_df_ix$treatment = as.factor(all_treatments_df_ix$treatment);

  all_treatments_plot_ix = (
    all_treatments_df_ix %>%
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
    scale_color_manual(values=line_colors) +
    scale_fill_manual(values=line_colors) +
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
    all_treatments_plot_ix,
    width=18,
    height=9);
}
