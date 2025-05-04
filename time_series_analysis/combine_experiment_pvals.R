library(dplyr)
library(reshape2)


cgrpaba_power_pval_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/CGRPaba",
  "/power",
  "/gam_analysis",
  "/CGRPaba_GAM_analysis_pvalues.csv");

cgrpaba_coherence_pval_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/CGRPaba",
  "/coherence",
  "/gam_analysis",
  "/CGRPaba_GAM_analysis_pvalues.csv");

sumatriptan_power_pval_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/sumatriptan",
  "/power",
  "/gam_analysis",
  "/sumatriptan_GAM_analysis_pvalues.csv");

sumatriptan_coherence_pval_filename = paste0(
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

els_power_pval_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/ELS",
  "/power",
  "/gam_analysis",
  "/ELS_GAM_analysis_pvalues.csv");

els_coherence_pval_filename = paste0(
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

combined_pval_save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/brainbo_lab_repo",
  "/brainbolab-code-base",
  "/time_series_analysis",
  "/combined_pvals_across_experiments_w_fdr_adjustments.csv");

significance_level = 0.05;

pow_melt_id_vars = c("region", "freq_band");
coh_melt_id_vars = c("region_pair", "freq_band");

combined_pval_df = (
  rbind(
    read.csv(cgrpaba_power_pval_filename, header=TRUE) %>%
    melt(id.vars=pow_melt_id_vars) %>%
    rename(pval_type=variable, pval=value) %>%
    mutate(experiment="CGRPaba", data_type="power"),

    read.csv(cgrpaba_coherence_pval_filename, header=TRUE) %>%
    melt(id.vars=coh_melt_id_vars) %>%
    rename(region=region_pair, pval_type=variable, pval=value) %>%
    mutate(experiment="CGRPaba", data_type="coherence"),

    read.csv(sumatriptan_power_pval_filename, header=TRUE) %>%
    melt(id.vars=pow_melt_id_vars) %>%
    rename(pval_type=variable, pval=value) %>%
    mutate(experiment="sumatriptan", data_type="power"),

    read.csv(sumatriptan_coherence_pval_filename, header=TRUE) %>%
    melt(id.vars=coh_melt_id_vars) %>%
    rename(region=region_pair, pval_type=variable, pval=value) %>%
    mutate(experiment="sumatriptan", data_type="coherence"),

    read.csv(els_power_pval_filename, header=TRUE) %>%
    melt(id.vars=pow_melt_id_vars) %>%
    rename(pval_type=variable, pval=value) %>%
    mutate(experiment="ELS", data_type="power"),

    read.csv(els_coherence_pval_filename, header=TRUE) %>%
    melt(id.vars=coh_melt_id_vars) %>%
    rename(region=region_pair, pval_type=variable, pval=value) %>%
    mutate(experiment="ELS", data_type="coherence") ) %>%

  mutate(
    pval_fdr_adjusted=p.adjust(pval, method="fdr"),
    fdr_significant=pval_fdr_adjusted < significance_level) );

order_ixs = order(
  combined_pval_df$experiment,
  combined_pval_df$data_type,
  combined_pval_df$region,
  combined_pval_df$freq_band);

combined_pval_df = combined_pval_df[order_ixs,];
rownames(combined_pval_df) = NULL;

write.csv(combined_pval_df, combined_pval_save_filename, row.names=FALSE);
