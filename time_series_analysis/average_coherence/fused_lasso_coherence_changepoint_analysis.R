library(genlasso)
library(ggplot2)
library(parallel)
library(purrr)
library(stringr)


###################################################################################
## Change this to the path of the directory where the mean coherence time series ##
## table is saved.                                                               ##
###################################################################################
mean_coherence_ts_filename = paste0(
  "/Users",
  "/ikhultman",
  "/power_time_series_analysis",
  "/CGRP_power_time_series_analysis",
  "/fused_lasso_analysis",
  "/coherence_changepoint_analysis",
  "/fused_lasso_coherence_data240408",
  "/Day3_coherence_fused_lasso_analysis.csv");

#################################################################################
## Change this to the path of the directory where the results of this analysis ##
## should be saved.                                                            ##
#################################################################################
save_dir = paste0(
  "/Users",
  "/ikhultman",
  "/power_time_series_analysis",
  "/CGRP_power_time_series_analysis",
  "/fused_lasso_analysis",
  "/coherence_changepoint_analysis",
  "/fused_lasso_coherence_results240408",
  "/Day3_results",
  "/post_injection/");

if (!dir.exists(save_dir) ) {
  dir.create(save_dir, recursive=TRUE);

  stopifnot(dir.exists(save_dir) );
}

# Location where the plots constructed by this fused lasso analysis will be saved.
figures_dir = paste0(save_dir, "figures/");

if (!dir.exists(figures_dir) ) {
  dir.create(figures_dir);

  stopifnot(dir.exists(figures_dir) );
}

# Location where the plots showing how mean test MSE changes with increasing
# lambda values will be saved. These plots also show a one standard error region
# around the smallest mean test MSE which is used to find the largest lambda
# value that provides decent prediction results.
figures_lam_choice_dir = paste0(figures_dir, "lambda_choice/");

if (!dir.exists(figures_lam_choice_dir) ) {
  dir.create(figures_lam_choice_dir);

  stopifnot(dir.exists(figures_lam_choice_dir) );
}

# Location where the plots comparing the observed data to the fused lasso
# fitted values for each brain region will be saved.
figures_fit_v_obs_dir = paste0(figures_dir, "fitted_vs_observed/");

if (!dir.exists(figures_fit_v_obs_dir) ) {
  dir.create(figures_fit_v_obs_dir);

  stopifnot(dir.exists(figures_fit_v_obs_dir) );
}

lambda_1se_dir = paste0(figures_fit_v_obs_dir, "lambda_1se");

if (!dir.exists(lambda_1se_dir) ) {
  dir.create(lambda_1se_dir);
  stopifnot(dir.exists(lambda_1se_dir) );
}

lambda_min_dir = paste0(figures_fit_v_obs_dir, "lambda_min");

if (!dir.exists(lambda_min_dir) ) {
  dir.create(lambda_min_dir);
  stopifnot(dir.exists(lambda_min_dir) );
}

mean_coherence_ts_data = read.csv(mean_coherence_ts_filename, header=TRUE);

expected_noncoh_colnames = c("mouse_id", "freq_band", "region_pair");
noncoh_col_ixs = which(
  is.element(colnames(mean_coherence_ts_data), expected_noncoh_colnames) );

stopifnot(length(noncoh_col_ixs) == length(expected_noncoh_colnames) );

noncoher_info = mean_coherence_ts_data[,noncoh_col_ixs];
n_obs = nrow(noncoher_info);

#segment_range_ixs = c(1, 1801);
segment_range_ixs = c(1802, 5403);

coherence_col_ixs = sort(
  setdiff(
    1:ncol(mean_coherence_ts_data),
    noncoh_col_ixs) )[segment_range_ixs[1]:segment_range_ixs[2]];

mean_coherence_ts_data = mean_coherence_ts_data[,c(noncoh_col_ixs, coherence_col_ixs)];
coherence_col_ixs = (length(noncoh_col_ixs) + 1):ncol(mean_coherence_ts_data);

coherence_mat_init = as.matrix(mean_coherence_ts_data[,coherence_col_ixs]);
thinning_factor = 15;

n_blocks = ceiling(length(coherence_col_ixs) / thinning_factor);
block_sizes = rep(thinning_factor, n_blocks);
n_remaining = length(coherence_col_ixs) %% thinning_factor;

if (n_remaining > 0) {
  block_sizes[n_blocks] = n_remaining;
}

stopifnot(sum(block_sizes) == length(coherence_col_ixs) );

block_ixs_table = rbind(
  cumsum(c(1, block_sizes[1:(n_blocks-1)]) ),
  cumsum(block_sizes) );

coherence_mat = matrix(NA, n_obs, n_blocks);

for (bx in 1:n_blocks) {
  col_ax = block_ixs_table[1,bx];
  col_bx = block_ixs_table[2,bx];

  if (block_sizes[bx] > 1) {
    coherence_mat[,bx] = rowMeans(coherence_mat_init[,col_ax:col_bx], na.rm=TRUE);
  } else {
    coherence_mat[,bx] = coherence_mat_init[,col_ax];
  }
}

mean_coherence_ts_data = cbind(noncoher_info, coherence_mat);
coherence_col_ixs = (1:n_blocks) + ncol(noncoher_info);

#coherence_mat_test = coherence_mat_init;

#if (n_remaining > 0) {
#  coherence_mat_test = cbind(
#    coherence_mat_test,
#    matrix(NA, n_obs, thinning_factor - n_remaining) );
#}

#n_coherence_pts = ncol(coherence_mat_test) / thinning_factor;
#coherence_arr = array(t(coherence_mat_test), dim=c(thinning_factor, n_coherence_pts, n_obs) );
#coherence_mat_test = apply(
#  coherence_arr,
#  MARGIN=c(2, 3),
#  FUN=function(nxt_block) mean(nxt_block, na.rm=TRUE) );

#stopifnot(mean((t(coherence_mat_test) - coherence_mat)^2) == 0);

#coherence_col_ixs = seq(
#  coherence_col_ixs[1],
#  coherence_col_ixs[length(coherence_col_ixs)],
#  thinning_factor);

# Construct the D penalty matrix for fused lasso regression.
D = getDtf(n_blocks, 0);

# Plots background color.
bg_color = rgb(245, 245, 245, maxColorValue=255);

n_regpair_freq_combos = nrow(unique(mean_coherence_ts_data[,c("region_pair", "freq_band")]) );

# The mean fitted values at each frequency for each pair of brain regions will
# be stored in the following tables and saved to a CSV file.
results_mean_fitted_1se = as.data.frame(matrix(NA, n_blocks, n_regpair_freq_combos) );
results_mean_fitted_min = as.data.frame(matrix(NA, n_blocks, n_regpair_freq_combos) );

# The standard deviation of fitted values at each frequency for each pair of brain
# regions will be stored in the following tables and saved to a CSV file.
results_sd_fitted_1se = as.data.frame(matrix(NA, n_blocks, n_regpair_freq_combos) );
results_sd_fitted_min = as.data.frame(matrix(NA, n_blocks, n_regpair_freq_combos) );

# The mean observed values at each frequency for each pair of brain regions will
# be stored in the following table and saved to a CSV file.
results_mean_observed = as.data.frame(matrix(NA, n_blocks, n_regpair_freq_combos) );

brain_region_pairs = unique(mean_coherence_ts_data$region_pair);
n_region_pairs = length(brain_region_pairs);

col_ix = 0;


fused_lasso_loocv = function(Y, X, D, log_lambdas, gamma) {
# DESCRIPTION:
#   Leave-one-out cross-validation for fused lasso.
#
# PARAMETERS:
#             Y: Response vector.
#
#             X: Design matrix for fused lasso.
#
#             D: Fused lasso penalty matrix.
#
#   log_lambdas: If log-lambda(s) is(are) provided, then that will be used for
#                cross-validation tests; otherwise log-lambdas over the range from
#                log(1e-6) to the largest produced during training are used.
#
#         gamma: Parameter used to penalize beta parameter estimates such that
#                larger gamma values favor sparser models.
#
# RETURN VALUE:
#   LOOCV test MSE means, standard deviations and corresponding log-lambda values.

  p_dim = ncol(X);
  n_dim = nrow(X) / p_dim;

  tst_inds = matrix(1:(p_dim*n_dim), ncol=n_dim);
  trn_inds = apply(
    tst_inds,
    MARGIN=2,
    FUN=function(col) setdiff(1:(p_dim*n_dim), col) );

  cat(sprintf("Fitting %d leave-one-out cross-validation models ...\n", n_dim) );

  mods_train = lapply(
    as.list(1:n_dim),
    FUN=function(gx) {
      X_train = X[trn_inds[,gx],];
      Y_train = Y[trn_inds[,gx]];
      fusedlasso(Y_train, X_train, D, gamma=gamma)
    });

  if (any(is.na(log_lambdas) )) {
    n_lambdas = 100;
    log_lam_min = log(1e-6);
    log_lam_max = log(
      ceiling(
        max(
          unlist(
            lapply(
              mods_train,
              FUN=function(mod) max(mod$lambda) )))));

    log_lambdas = seq(log_lam_min, log_lam_max, length.out=n_lambdas);
  }

  n_lambdas = length(log_lambdas);

  cat(
    sprintf(
      "Computing LOOCV test MSE's at %d different lambda values ...\n\n",
      n_lambdas) );

  loocv_results = data.frame(
    log_lambda=log_lambdas,
    mean_mse=rep(NA, n_lambdas),
    sd_mse=rep(NA, n_lambdas) );

  all_mse = matrix(NA, n_lambdas, n_dim);

  for (gx in 1:n_dim) {
    next_mod_trn = mods_train[[gx]];
    betas = coef(next_mod_trn, lambda=exp(loocv_results$log_lambda) )$beta;
    X_test = X[tst_inds[,gx],];
    Y_test = Y[tst_inds[,gx]];
    Y_hat = X_test %*% betas;
    all_mse[,gx] = colMeans((Y_test - Y_hat)^2);
  }

  loocv_results$mean_mse = rowMeans(all_mse);
  loocv_results$sd_mse = apply(all_mse, 1, sd);

  list(
    loocv_results=loocv_results,
    mods_train=mods_train);
};


fused_lasso_loocv_par = function(Y, X, D, log_lambdas, gamma, n_proc) {
# DESCRIPTION:
#   Leave-one-out cross-validation for fused lasso.
#
# PARAMETERS:
#             Y: Response vector.
#
#             X: Design matrix for fused lasso.
#
#             D: Fused lasso penalty matrix.
#
#   log_lambdas: If log-lambda(s) is(are) provided, then that will be used for
#                cross-validation tests; otherwise log-lambdas over the range from
#                log(1e-6) to the largest produced during training are used.
#
#         gamma: Parameter used to penalize beta parameter estimates such that
#                larger gamma values favor sparser models.
#
#            cl: Cluster object for parallel computation.
#
# RETURN VALUE:
#   LOOCV test MSE means, standard deviations and corresponding log-lambda values.

  p_dim = ncol(X);
  n_dim = nrow(X) / p_dim;

  tst_inds = matrix(1:(p_dim*n_dim), ncol=n_dim);
  trn_inds = apply(
    tst_inds,
    MARGIN=2,
    FUN=function(col) setdiff(1:(p_dim*n_dim), col) );

  cat(sprintf("Fitting %d leave-one-out cross-validation models ...\n", n_dim) );

  cl_obj = makeCluster(min(n_proc, n_dim) );
  cl_lib_info = clusterEvalQ(cl_obj, library(genlasso) );
  clusterExport(cl_obj, c("X", "Y", "D", "gamma", "trn_inds"), environment() );

  mods_train = parLapply(
    cl_obj,
    as.list(1:n_dim),
    function(gx) {
      X_train = X[trn_inds[,gx],];
      Y_train = Y[trn_inds[,gx]];
      fusedlasso(Y_train, X_train, D, gamma=gamma)
    });

  stopCluster(cl_obj);

  if (any(is.na(log_lambdas) )) {
    n_lambdas = 100;
    log_lam_min = log(1e-6);
    log_lam_max = log(
      ceiling(
        max(
          unlist(
            lapply(
              mods_train,
              FUN=function(mod) max(mod$lambda) )))));

    log_lambdas = seq(log_lam_min, log_lam_max, length.out=n_lambdas);
  }

  n_lambdas = length(log_lambdas);

  cat(
    sprintf(
      "Computing LOOCV test MSE's at %d different lambda values ...\n\n",
      n_lambdas) );

  loocv_results = data.frame(
    log_lambda=log_lambdas,
    mean_mse=rep(NA, n_lambdas),
    sd_mse=rep(NA, n_lambdas) );

  all_mse = matrix(NA, n_lambdas, n_dim);

  for (gx in 1:n_dim) {
    next_mod_trn = mods_train[[gx]];
    betas = coef(next_mod_trn, lambda=exp(loocv_results$log_lambda) )$beta;
    X_test = X[tst_inds[,gx],];
    Y_test = Y[tst_inds[,gx]];
    Y_hat = X_test %*% betas;
    all_mse[,gx] = colMeans((Y_test - Y_hat)^2);
  }

  loocv_results$mean_mse = rowMeans(all_mse);
  loocv_results$sd_mse = apply(all_mse, 1, sd);

  list(
    loocv_results=loocv_results,
    mods_train=mods_train);
};


for (rx in 1:n_region_pairs) {
  next_region_pair = brain_region_pairs[rx];
  next_region_pair_ixs = mean_coherence_ts_data$region_pair == next_region_pair;
  next_region_pair_data = mean_coherence_ts_data[next_region_pair_ixs,];
  rownames(next_region_pair_data) = NULL;

  freq_bands = unique(next_region_pair_data$freq_band);
  n_freq_bands = length(freq_bands);

  for (fx in 1:n_freq_bands) {
    col_ix = col_ix + 1;

    next_freq_band = freq_bands[fx];
    next_freq_ixs = next_region_pair_data$freq_band == next_freq_band;
    next_reg_freq_data = next_region_pair_data[next_freq_ixs,];
    rownames(next_reg_freq_data) = NULL;

    n_mice = length(unique(next_reg_freq_data$mouse_id) );

    stopifnot(nrow(next_reg_freq_data) == n_mice);

    # Construct the X design matrix for fused lasso regression.
    X = diag(n_blocks);
    X = t(matrix(rep(X, n_mice), n_blocks) );

    #y_scaled_mat = scale(t(as.matrix(next_reg_freq_data[,coherence_col_ixs]) ));
    #Y = as.vector(y_scaled_mat);

    Y = as.vector(t(as.matrix(next_reg_freq_data[,coherence_col_ixs]) ));

    gamma = 0;
    loocv_results = fused_lasso_loocv_par(Y, X, D, NA, gamma, n_mice);

    mods_train = loocv_results$mods_train;
    loocv_results = loocv_results$loocv_results;

    min_tst_mse = min(loocv_results$mean_mse);
    min_mse_ix = max(which(loocv_results$mean_mse == min_tst_mse) );
    tst_mse_1se = min_tst_mse + ((loocv_results$sd_mse[min_mse_ix] / sqrt(n_mice) ) * c(-1, 1) );
    lam_1se_ix = max(which(loocv_results$mean_mse <= tst_mse_1se[2]) );
    log_lam_1se = loocv_results$log_lambda[lam_1se_ix];
    log_lam_min = loocv_results$log_lambda[min_mse_ix];

    plot_df1 = data.frame(
      log_lambda=loocv_results$log_lambda,
      mean_mse=loocv_results$mean_mse,
      min_tst_mse=rep(min_tst_mse, nrow(loocv_results) ),
      mse_1se_lower=rep(tst_mse_1se[1], nrow(loocv_results) ),
      mse_1se_upper=rep(tst_mse_1se[2], nrow(loocv_results) ));

    rm_ixs = is.infinite(plot_df1$log_lambda) | is.na(plot_df1$log_lambda);
    plot_df1 = plot_df1[!rm_ixs,];

    plot_df2 = data.frame(
      log_lambda_min=loocv_results$log_lambda[min_mse_ix],
      log_lambda_1se=loocv_results$log_lambda[lam_1se_ix]);

    next_mse_results_plot = (
      ggplot() +
      theme(
        panel.background=element_rect(fill=bg_color) ) +
      geom_line(
        data=plot_df1,
        aes(
          x=log_lambda,
          y=mean_mse,
          colour="Mean Test MSE"),
        linewidth=1) +
      geom_vline(
        data=plot_df2,
        aes(
          xintercept=log_lambda_min,
          colour="Log(lambda_min)"),
        linetype="dotted",
        linewidth=0.7,
        show.legend=FALSE) +
      geom_vline(
        data=plot_df2,
        aes(
          xintercept=log_lambda_1se,
          colour="Log(lambda_1se)"),
        linetype="dashed",
        linewidth=0.7,
        show.legend=FALSE) +
      geom_ribbon(
        data=plot_df1,
        aes(
          x=log_lambda,
          y=min_tst_mse,
          ymin=mse_1se_lower,
          ymax=mse_1se_upper,
          colour="1 Std. Error"),
        alpha=0.2,
        show.legend=FALSE) +
      scale_colour_manual(
        values=c("gray", "red", "blue", "black"),
        guide=guide_legend(
          title="Legend",
          override.aes=list(
            fill=c("gray40", bg_color, bg_color, bg_color),
            linetype=c("solid", "dashed", "dotted", "solid"),
            linewidth=c(6, 0.7, 0.7, 1) ))) +
      labs(
        x="Log(lambda)",
        y="Mean Test MSE") );

    save_filename = paste0(
      figures_lam_choice_dir,
      sprintf(
        "lambda_choice_coherence_%s_freq_band_%d.png",
        next_region_pair,
        next_freq_band) );

    ggsave(
      save_filename,
      next_mse_results_plot,
      units="in",
      height=7,
      width=7,
      device="png");

    #mod_full_data = fusedlasso(Y, X, D, gamma=gamma);

    results_colname = paste0(next_region_pair, "_freq_band_", next_freq_band);

    mods_train_fitted_1se = matrix(
      unlist(
        lapply(
          mods_train,
          function(next_mod_trn) as.numeric(coef(next_mod_trn, lambda=exp(log_lam_1se) )$beta) )),
      ncol=n_mice);

    results_mean_fitted_1se[,col_ix] = rowMeans(mods_train_fitted_1se);
    colnames(results_mean_fitted_1se)[col_ix] = results_colname;

    results_sd_fitted_1se[,col_ix] = apply(mods_train_fitted_1se, 1, sd);
    colnames(results_sd_fitted_1se)[col_ix] = results_colname;

    mods_train_fitted_min = matrix(
      unlist(
        lapply(
          mods_train,
          function(next_mod_trn) as.numeric(coef(next_mod_trn, lambda=exp(log_lam_min) )$beta) )),
      ncol=n_mice);

    results_mean_fitted_min[,col_ix] = rowMeans(mods_train_fitted_min);
    colnames(results_mean_fitted_min)[col_ix] = results_colname;

    results_sd_fitted_min[,col_ix] = apply(mods_train_fitted_min, 1, sd);
    colnames(results_sd_fitted_min)[col_ix] = results_colname;

    #results_fitted_1se[,col_ix] = as.numeric(coef(mod_full_data, lambda=exp(log_lam_1se) )$beta);
    #colnames(results_fitted_1se)[col_ix] = paste0(next_region_pair, "_freq_band_", next_freq_band);

    #results_fitted_min[,col_ix] = as.numeric(coef(mod_full_data, lambda=exp(log_lam_min) )$beta);
    #colnames(results_fitted_min)[col_ix] = paste0(next_region_pair, "_freq_band_", next_freq_band);

    results_mean_observed[,col_ix] = as.numeric(colMeans(as.matrix(next_reg_freq_data[,coherence_col_ixs]) ));
    colnames(results_mean_observed)[col_ix] = results_colname;

    save_filename = sprintf(
      "fused_lasso_lambda_1se_coherence_%s_freq_band_%d.png",
      next_region_pair,
      next_freq_band);

    save_filename = paste0(
      figures_fit_v_obs_dir,
      "lambda_1se/",
      save_filename);

    cat(
      sprintf(
        paste0(
          "Saving plot comparing fitted and observed values for ",
          "lambda w/ 1 std. err. test MSE to:\n%s\n\n"),
        save_filename) );

    png(filename=save_filename);

    plot_xs = ((0:(n_blocks - 1) ) * thinning_factor) + 1;

    #mean_scaled_y = as.numeric(rowMeans(y_scaled_mat) );

    mean_y = results_mean_observed[,col_ix];

    #y_lim = c(
    #  min(-1, min(mean_scaled_y) ),
    #  max(1, max(mean_scaled_y) ));

    y_lim = c(min(mean_y), max(mean_y) );
    y_lim_range = diff(y_lim);
    y_scale_inc = 0.2;
    y_inc = y_lim_range * y_scale_inc;
    y_lim = y_lim + ((0.5 * y_inc) * c(-1, 1) );

    plot(
      plot_xs,
      mean_y,
      lty=2,
      col="black",
      main=paste0(
        next_region_pair,
        " - frequency band ",
        next_freq_band),
      ylab="Coherence",
      xlab="Time (seconds)",
      type="l",
      ylim=y_lim);

    lines(plot_xs, results_mean_fitted_1se[[col_ix]], col="steelblue", lwd=3);

    #abline(v=1801, col="green4", lty=6, lwd=2);
    #legend(
    #  0.65 * n_blocks * thinning_factor, 1,
    #  legend=c("Mean Observed", "Fitted", "Injection"),
    #  col=c("black", "steelblue", "green4"),
    #  lty=c(2, 1, 6),
    #  lwd=c(1, 3, 2) );

    legend(
      0.65 * n_blocks * thinning_factor, y_lim[2],
      legend=c("Mean Observed", "Fitted"),
      col=c("black", "steelblue"),
      lty=c(2, 1),
      lwd=c(1, 3) );

    dev.off();

    save_filename = sprintf(
      "fused_lasso_lambda_min_coherence_%s_freq_band_%d.png",
      next_region_pair,
      next_freq_band);

    save_filename = paste0(
      figures_fit_v_obs_dir,
      "lambda_min/",
      save_filename);

    cat(
      sprintf(
        paste0(
          "Saving plot comparing fitted and observed values for ",
          "lambda w/ min test MSE to:\n%s\n\n"),
        save_filename) );

    png(filename=save_filename);

    plot(
      plot_xs,
      mean_y,
      lty=2,
      col="black",
      main=paste0(
        next_region_pair,
        " - frequency band ",
        next_freq_band),
      ylab="Coherence",
      xlab="Time (seconds)",
      type="l",
      ylim=y_lim);

    lines(plot_xs, results_mean_fitted_min[[col_ix]], col="steelblue", lwd=3);

    #abline(v=1801, col="green4", lty=6, lwd=2);
    #legend(
    #  0.65 * n_blocks * thinning_factor, 1,
    #  legend=c("Mean Observed", "Fitted", "Injection"),
    #  col=c("black", "steelblue", "green4"),
    #  lty=c(2, 1, 6),
    #  lwd=c(1, 3, 2) );

    legend(
      0.65 * n_blocks * thinning_factor, y_lim[2],
      legend=c("Mean Observed", "Fitted"),
      col=c("black", "steelblue"),
      lty=c(2, 1),
      lwd=c(1, 3) );

    dev.off();
  }
}

save_filename = paste0(
  save_dir,
  "mean_fitted_values_lambda_1se.csv");

write.table(
  results_mean_fitted_1se,
  save_filename,
  sep=", ",
  row.names=FALSE);

save_filename = paste0(
  save_dir,
  "sd_fitted_values_lambda_1se.csv");

write.table(
  results_sd_fitted_1se,
  save_filename,
  sep=", ",
  row.names=FALSE);

save_filename = paste0(
  save_dir,
  "mean_fitted_values_lambda_min.csv");

write.table(
  results_mean_fitted_min,
  save_filename,
  sep=", ",
  row.names=FALSE);

save_filename = paste0(
  save_dir,
  "sd_fitted_values_lambda_min.csv");

write.table(
  results_sd_fitted_min,
  save_filename,
  sep=", ",
  row.names=FALSE);

save_filename = paste0(
  save_dir,
  "mean_observed_values.csv");

write.table(
  results_mean_observed,
  save_filename,
  sep=", ",
  row.names=FALSE);
