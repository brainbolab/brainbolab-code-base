# This script takes the path to a CSV file that contains a table of the following form:
#
#      ID   date Treatment Sex phase1_start phase1_end phase2_start phase2_end ...
#   MV211 240517   CGRP_A1   M            1        600          862       1019 ...
#   MV212 240517   CGRP_A1   M            1        600          847       1004 ...
#   MV214 240605   CGRP_A1   M            1        600          832        989 ...
#   MV217 240613   CGRP_A1   M            1        600          849       1006 ...
#   MV218 240621   CGRP_A1   M            1        600          815        972 ...
#   MV221 240621   CGRP_A1   F            1        600          831        988 ...
#   ...      ...       ... ...          ...        ...          ...        ... ...
#
# and stretches it according to the experimental phases denoted by columns that end in '_start'
# and '_end'. This results in a new table that has the following form:
#
#      ID   date Treatment Sex exp_phase start_time  end_time ...
#   MV211 240517   CGRP_A1   M    phase1          1       600 ...
#   MV211 240517   CGRP_A1   M    phase2        862      1019 ...
#     ...    ...       ... ...       ...        ...       ... ...
#   MV212 240517   CGRP_A1   M    phase1          1       600 ...
#   MV212 240517   CGRP_A1   M    phase2        847      1004 ...
#     ...    ...       ... ...       ...        ...       ... ...
#   MV214 240605   CGRP_A1   M    phase1          1       600 ...
#   MV214 240605   CGRP_A1   M    phase2        832       989 ...
#     ...    ...       ... ...       ...        ...       ... ...
#   MV217 240613   CGRP_A1   M    phase1          1       600 ...
#   MV217 240613   CGRP_A1   M    phase2        849      1006 ...
#     ...    ...       ... ...       ...        ...       ... ...
#   MV218 240621   CGRP_A1   M    phase1          1       600 ...
#   MV218 240621   CGRP_A1   M    phase2        815       972 ...
#     ...    ...       ... ...       ...        ...       ... ...
#   MV221 240621   CGRP_A1   F    phase1          1       600 ...
#   MV221 240621   CGRP_A1   F    phase2        831       988 ...
#     ...    ...       ... ...       ...        ...       ... ...
#
# to which we can further label each time segment.

library(stringr)
library(tidyr)


# Path to the original times table.
times_table_fname = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/light_dark_cvsfa",
  "/light_dark_CVSFA",
  "/lpne_analysis",
  "/yassine_experiments",
  "/MigraineV2xFLA_C2_C3_times_tables_CGRP_A1_A2_B_250324.csv");

# Path where the stretched times and labels table will be saved.
save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/light_dark_cvsfa",
  "/light_dark_CVSFA",
  "/lpne_analysis",
  "/yassine_experiments",
  "/MigraineV2xFLA_CGRP_times_and_labels.csv");

times_table = read.csv(times_table_fname, header=TRUE);
n_rec = nrow(times_table);

exp_phases_w_start = str_extract(colnames(times_table), ".*(?=_start$)");
exp_phases_w_start = exp_phases_w_start[!is.na(exp_phases_w_start)];

exp_phases_w_end = str_extract(colnames(times_table), ".*(?=_end$)");
exp_phases_w_end = exp_phases_w_end[!is.na(exp_phases_w_end)];

exp_phases = intersect(exp_phases_w_start, exp_phases_w_end);
n_exp_phases = length(exp_phases);

time_labels = as.vector(
  sapply(
    exp_phases,
    \(phase_ix) paste0(phase_ix, c("_start", "_end") )));

non_time_labels = setdiff(colnames(times_table), time_labels);

times_table$time_regions = vector(mode="list", length=n_rec);

for (rx in 1:n_rec) {
  times_table$time_regions[[rx]] = data.frame(
    exp_phase=character(n_exp_phases),
    start_time=NA,
    end_time=NA);

  for (px in 1:n_exp_phases) {
    phase_px = exp_phases[px];
    times_table$time_regions[[rx]][px,"exp_phase"] = phase_px;
    times_table$time_regions[[rx]][px,"start_time"] = times_table[rx,paste0(phase_px, "_start")];
    times_table$time_regions[[rx]][px,"end_time"] = times_table[rx,paste0(phase_px, "_end")];
  }
}

times_table = unnest(times_table[,c(non_time_labels, "time_regions")], cols="time_regions");

# Adding CGRP vs vehicle labels to the time segments that are either CGRP or vehicle.
times_table$CGRP = "veh";
ixs_cgrp = (times_table$exp_phase != "baseline") & (times_table$Treatment == "CGRP_B");
times_table$CGRP[ixs_cgrp] = "CGRP";

write.csv(times_table, save_filename, row.names=FALSE);
