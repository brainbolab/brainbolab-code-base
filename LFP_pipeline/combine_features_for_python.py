import numpy as np
import os
import pandas as pd
import pickle
import re
from scipy.io import loadmat
import warnings


data_dir = "./LFP_features";
save_filename = "./LFP_py_features.pickle";

if not os.path.isdir(data_dir):
  err_msg = "The provided directory:\n  " + data_dir + "\ndoes not exist.";
  raise AssertionError(err_msg);

filenames = sorted(os.listdir(data_dir) );
n_filenames = len(filenames);

re_mouse_id = "(?i)((?<=_mouse)[a-z0-9]+(?=_|\\.))|((?<=^mouse)[a-z0-9]+(?=_|\\.))";
re_date = "((?<=_)[0-9]{6}(?=_|\\.))|((?<=^)[0-9]{6}(?=_|\\.))";

combined_features = {};

for fx in range(0, n_filenames):
  print("Checking features for file " + filenames[fx] + " ...");

  if os.path.splitext(filenames[fx])[1] != ".mat":
    warn_msg = "The file " + filenames[fx] + " does not have .mat extension. Skipping ...";
    warnings.warn(warn_msg);
    continue;

  fname_fx = os.path.join(data_dir, filenames[fx]);
  features_fx = loadmat(fname_fx);

  power_fx = features_fx["power"];
  fft_fx = features_fx["fft"];
  coherence_fx = features_fx["coherence"];
  labels_fx = [labl_ix[0][0] for labl_ix in features_fx["labels"]];
  regions_fx = np.array([reg_ix[0][0] for reg_ix in features_fx["regions"]]);
  region_pairs_fx = np.array([reg_ix[0][0] for reg_ix in features_fx["region_pairs"]]);
  freq_fx = features_fx["freq"].flatten();

  n_win_fx = len(labels_fx);
  n_regions_fx = len(regions_fx);
  n_reg_pairs_fx = len(region_pairs_fx);
  n_freq_fx = len(freq_fx);

  if len(power_fx.shape) != 3:
    warn_msg = "The power array from file " + filenames[fx] + " is not 3D. Skipping ...";
    warnings.warn(warn_msg);
    continue;

  if len(fft_fx.shape) != 3:
    warn_msg = "The FFT array from file " + filenames[fx] + " is not 3D. Skipping ...";
    warnings.warn(warn_msg);
    continue;

  if len(coherence_fx.shape) != 3:
    warn_msg = "The coherence array from file " + filenames[fx] + " is not 3D. Skipping ...";
    warnings.warn(warn_msg);
    continue;

  if n_freq_fx != power_fx.shape[0]:
    warn_msg = \
      "The number of frequencies represented by the power array in file " +\
      filenames[fx] +\
      " does not equal the number of frequency labels. Skipping ...";

    warnings.warn(warn_msg);
    continue;

  if n_freq_fx != fft_fx.shape[0]:
    warn_msg = \
      "The number of frequencies represented by the FFT array in file " +\
      filenames[fx] +\
      " does not equal the number of frequency labels. Skipping ...";

    warnings.warn(warn_msg);
    continue;

  if n_freq_fx != coherence_fx.shape[0]:
    warn_msg = \
      "The number of frequencies represented by the coherence array in file " +\
      filenames[fx] +\
      " does not equal the number of frequency labels. Skipping ...";

    warnings.warn(warn_msg);
    continue;

  if n_regions_fx != power_fx.shape[1]:
    warn_msg = \
      "The number of regions represented by the power array in file " +\
      filenames[fx] +\
      " does not equal the number of region labels. Skipping ...";

    warnings.warn(warn_msg);

    continue;

  if n_regions_fx != fft_fx.shape[1]:
    warn_msg = \
      "The number of regions represented by the FFT array in file " +\
      filenames[fx] +\
      " does not equal the number of region labels. Skipping ...";

    warnings.warn(warn_msg);

    continue;

  if n_reg_pairs_fx != coherence_fx.shape[1]:
    warn_msg = \
      "The number of region pairs represented by the coherence array in file " +\
      filenames[fx] +\
      " does not equal the number of region pair labels. Skipping ...";

    warnings.warn(warn_msg);

    continue;

  if n_win_fx != power_fx.shape[2]:
    warn_msg = \
      "The number of windows represented by the power array in file " +\
      filenames[fx] +\
      " does not equal the number of window labels. Skipping ...";

    warnings.warn(warn_msg);
    continue;

  if n_win_fx != fft_fx.shape[2]:
    warn_msg = \
      "The number of windows represented by the FFT array in file " +\
      filenames[fx] +\
      " does not equal the number of window labels. Skipping ...";

    warnings.warn(warn_msg);
    continue;

  if n_win_fx != coherence_fx.shape[2]:
    warn_msg = \
      "The number of windows represented by the coherence array in file " +\
      filenames[fx] +\
      " does not equal the number of window labels. Skipping ...";

    warnings.warn(warn_msg);
    continue;

  if n_regions_fx > 1:
    ixs_order = np.argsort(regions_fx);
    regions_fx = regions_fx[ixs_order];
    power_fx = power_fx[:,ixs_order,:];
    fft_fx = fft_fx[:,ixs_order,:];

  if n_reg_pairs_fx > 1:
    ixs_order = np.argsort(region_pairs_fx);
    region_pairs_fx = region_pairs_fx[ixs_order];
    coherence_fx = coherence_fx[:,ixs_order,:];

  mouse_fx = re.search(re_mouse_id, filenames[fx]).group();
  date_fx = re.search(re_date, filenames[fx]).group();
  group_info_fx = pd.DataFrame({
    "ID": [mouse_fx] * n_win_fx,
    "date": [date_fx] * n_win_fx,
    "label": labels_fx,
    "time": list(range(0, n_win_fx) )});

  if len(combined_features) == 0:
    combined_features["power"] = power_fx;
    combined_features["fft"] = fft_fx;
    combined_features["coherence"] = coherence_fx;
    combined_features["group_info"] = group_info_fx;
    combined_features["regions"] = regions_fx;
    combined_features["region_pairs"] = region_pairs_fx;
    combined_features["freq"] = freq_fx;
  else:
    regions_diff_new = set(regions_fx).difference(combined_features["regions"]);

    if len(regions_diff_new) > 0:
      warn_msg = "File " + filenames[fx] + " contains more regions than previous files. Skipping ...";
      warnings.warn(warn_msg);
      continue;

    regions_diff_prev = set(combined_features["regions"]).difference(regions_fx);

    if len(regions_diff_prev) > 0:
      warn_msg = "File " + filenames[fx] + " is missing regions contained in the previous files. Skipping ...";
      warnings.warn(warn_msg);
      continue;

    reg_pairs_diff_new = set(region_pairs_fx).difference(combined_features["region_pairs"]);

    if len(reg_pairs_diff_new) > 0:
      warn_msg = "File " + filenames[fx] + " contains more region pairs than previous files. Skipping ...";
      warnings.warn(warn_msg);
      continue;

    reg_pairs_diff_prev = set(combined_features["region_pairs"]).difference(region_pairs_fx);

    if len(reg_pairs_diff_prev) > 0:
      warn_msg = "File " + filenames[fx] + " is missing region pairs contained in the previous files. Skipping ...";
      warnings.warn(warn_msg);
      continue;

    if (freq_fx != combined_features["freq"]).any():
      warn_msg = "The frequencies in file " + filenames[fx] + " do not match the frequencies from previous files. Skipping ...";
      warnings.warn(warn_msg);
      continue;

    # If the program reaches this point in a given loop iteration, the following must be true.
    assert (regions_fx == combined_features["regions"]).all();
    assert (region_pairs_fx == combined_features["region_pairs"]).all();

    combined_features["power"] = np.concatenate(
      (combined_features["power"], power_fx),
      axis=2);

    combined_features["fft"] = np.concatenate(
      (combined_features["fft"], fft_fx),
      axis=2);

    combined_features["coherence"] = np.concatenate(
      (combined_features["coherence"], coherence_fx),
      axis=2);

    combined_features["group_info"] = pd.concat(
      (combined_features["group_info"], group_info_fx),
      axis=0).reset_index(drop=True);

with open(save_filename, "wb") as save_file:
  pickle.dump(combined_features, save_file, protocol=pickle.HIGHEST_PROTOCOL);
