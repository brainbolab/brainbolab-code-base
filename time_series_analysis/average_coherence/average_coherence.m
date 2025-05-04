% For each pair of brain regions from the set of provided brain regions, this script will compute
% the mean coherence time series' over the indicated frequency ranges for each mouse on each date
% provided. The results are saved as a MATLAB struct with fields identifying tables containing mean
% coherence time series in one Hz increments for different phases of an experiment.

% Change this to the path of the directory containing the mice coherence data directories.
data_dir = [ ...
  '/Users', ...
  '/ikhultman', ...
  '/average_spectograms', ...
  '/test_data_coherence/'];

% Change this to the path where the results should be saved.
save_dir = [ ...
  '/Users', ...
  '/ikhultman', ...
  '/average_spectograms', ...
  '/mean_coherence_results/'];

% Change this to the path of the CSV file containing the mice IDs, dates and times at which
% the data will be trimmed.
times_table_filename = [ ...
  '/Users', ...
  '/ikhultman', ...
  '/average_spectograms', ...
  '/CGRPaba_Timepoints_230824.csv'];

% Change the values of this array such that the first element of each row is a mouse ID, and
% the second element of each row is the date corresponding to the data for that mouse for which
% the mean coherences time series' will be computed. The mouse IDs in this array must be a subset
% of the mouse IDs provided in the times table.
mice_w_dates = [ ...
  {'M15'}, {'220204'}; ...
  {'M19'}, {'220811'}];

% Change this array to include the subset of brain regions for which the mean coherence time series'
% will be computed for each pair of brain regions. Make sure that the brain region labels included
% as part of the data filenames have the same format as those listed here.
brain_regions = { ...
  'BLA', ...
  'CeA', ...
  'LPBN', ...
  'MD_thal'};

% Change this to include each frequency range in Hz over which to compute the mean coherence time
% series' for each brain region for each mouse.
freq_ranges = [ ...
  {4:7}, ...
  {11:18}];

% Change this to include each phase of the experiment for which the mean coherence time series' will
% be computed. Each experiment phase label that appears here must have corresponding "start" and "end"
% labels in the provided times table; e.g. if one of the phase labels is "Baseline", then there must be
% "Baseline_start" and "Baseline_end" fields in the times table identifying the first and last indices
% into each mouse's total dataset that correspond to that phase of the experiment.
experiment_phases = { ...
  'Baseline', ...
  'Post_injection'};

if iscolumn(brain_regions)
  brain_regions = brain_regions';
end

if iscolumn(freq_ranges)
  freq_ranges = freq_ranges';
end

if iscolumn(experiment_phases)
  experiment_phases = experiment_phases';
end

file_sep = '/';
if ispc
  file_sep = '\';
end

if data_dir(end) ~= file_sep
  data_dir = [data_dir, file_sep];
end

if save_dir(end) ~= file_sep
  save_dir = [save_dir, file_sep];
end

addpath(['.', file_sep, 'average_coherence_support']);

exp_phase_colnames = cellfun( ...
  @(exp_phase) [{[exp_phase, '_start']}, {[exp_phase, '_end']}], ...
  experiment_phases, ...
  'UniformOutput', false);

exp_phase_colnames = [exp_phase_colnames{:}];

expected_colnames = [ ...
  {'Animal_ID'}, ...
  {'Date'}, ...
  exp_phase_colnames];

table_opts = detectImportOptions(times_table_filename);
[~, ~, expected_ixs] = intersect(table_opts.VariableNames, expected_colnames);

if isrow(expected_ixs)
  expected_ixs = expected_ixs';
end

% Make sure that the provided table contains the variable names provided in the "expected_colnames"
% variable defined above. It's fine if the table contains more than just the variable names provided
% (and order doesn't matter), but it must contain at least those expected variable names.
assert( ...
  length(expected_ixs) == length(expected_colnames) && ...
  all(sort(expected_ixs) == (1:length(expected_colnames) )') );

date_ix = find(strcmp(table_opts.VariableNames, 'Date') );

% Check to make sure there aren't multiple date columns in the provided table.
assert(length(date_ix) == 1);

table_opts.VariableTypes{date_ix} = 'char';
times_table = readtable(times_table_filename, table_opts);

mice_w_dates_table = cell2table( ...
  mice_w_dates, ...
  'VariableNames', [{'Animal_ID'}, {'Date'}]);

n_mice_w_dates = size(mice_w_dates_table, 1);
mice_w_dates_table.Data_dir = cell(n_mice_w_dates, 1);

all_region_pairs = [];
regions_re_pattern = ['(', char(join(brain_regions, '|') ), ')'];
mice_filenames_table = [];

for ix = 1:n_mice_w_dates
  mice_w_dates_table.Data_dir{ix} = [ ...
    data_dir, ...
    'Mouse', ...
    mice_w_dates_table.Animal_ID{ix}, ...
    file_sep, ...
    mice_w_dates_table.Date{ix}, ...
    file_sep, ...
    'Coherence', ...
    file_sep];

  mouse_ix_fnames = ls_filenames_w_pattern( ...
    mice_w_dates_table.Data_dir{ix}, ...
    [regions_re_pattern, '.*\.mat$']);

  mouse_ix_reg_pairs_counts = cellfun( ...
    @(fname) numel(regexp(fname, regions_re_pattern) ), ...
    mouse_ix_fnames);

  mouse_ix_fnames = mouse_ix_fnames(mouse_ix_reg_pairs_counts == 2);
  n_mouse_ix_fnames = numel(mouse_ix_fnames);
  mouse_ix_reg_pairs = cellfun( ...
    @(fname) join(regexp(fname, regions_re_pattern, 'match'), '_x_'), ...
    mouse_ix_fnames);

  mouse_ix_fnames_table = [ ...
    repelem(mice_w_dates_table.Animal_ID(ix), n_mouse_ix_fnames, 1), ...
    repelem(mice_w_dates_table.Date(ix), n_mouse_ix_fnames, 1), ...
    mouse_ix_reg_pairs, ...
    mouse_ix_fnames];

  if isempty(mice_filenames_table)
    mice_filenames_table = mouse_ix_fnames_table;
  else
    mice_filenames_table = [mice_filenames_table; mouse_ix_fnames_table];
  end

  all_region_pairs = union(all_region_pairs, unique(mouse_ix_reg_pairs) );
end

if iscolumn(all_region_pairs)
  all_region_pairs = all_region_pairs';
end

mice_filenames_table = cell2table( ...
  mice_filenames_table, ...
  'VariableNames', [{'Animal_ID'}, {'Date'}, {'region_pair'}, {'filename'}]);

mice_all_info_table = join(mice_w_dates_table, times_table);

n_region_pairs = numel(all_region_pairs);
n_freq_ranges = numel(freq_ranges);
n_freqs_each = arrayfun(@(rx) length(freq_ranges{rx}), 1:n_freq_ranges);

all_freqs = [freq_ranges{:}];
n_tot_freqs = length(all_freqs);

mice_means_ts_colnames = [ {'mouse'},   {'date'}, {'freq_band'},    {'Hz'}, {'region_pair'}];
mice_means_ts_vartypes = [{'string'}, {'string'},     {'int64'}, {'int64'},      {'string'}];

assert(numel(mice_means_ts_colnames) == numel(mice_means_ts_vartypes) );

mice_means_ts = table( ...
  'Size', [n_tot_freqs * n_mice_w_dates * n_region_pairs, numel(mice_means_ts_colnames)], ...
  'VariableTypes', mice_means_ts_vartypes, ...
  'VariableNames', mice_means_ts_colnames);

frq_grp_ids = arrayfun( ...
  @(rx) repelem(rx, n_freqs_each(rx) ), ...
  1:n_freq_ranges, ...
  'UniformOutput', false);

frq_grp_ids = repmat(repmat([frq_grp_ids{:}], 1, n_region_pairs), 1, n_mice_w_dates);
mice_means_ts.freq_band = frq_grp_ids';
mice_means_ts.Hz = repmat(all_freqs', n_region_pairs * n_mice_w_dates, 1);
mice_means_ts.region_pair = repmat( ...
  repelem(all_region_pairs, 1, n_tot_freqs)', n_mice_w_dates, 1);

n_exp_phases = numel(experiment_phases);
n_phase_pts_vec = NaN(1, n_exp_phases);

exp_phases_data = struct();

for px = 1:n_exp_phases
  exp_phase = experiment_phases{px};
  exp_phase_start = [exp_phase, '_start'];
  exp_phase_end = [exp_phase, '_end'];
  n_exp_phase_pts = mice_all_info_table.(exp_phase_end) - mice_all_info_table.(exp_phase_start) + 1;

  assert(all(n_exp_phase_pts == n_exp_phase_pts(1) ));

  n_phase_pts_vec(px) = n_exp_phase_pts(1);
  exp_phases_data.(exp_phase) = NaN( ...
    n_phase_pts_vec(px), ...
    n_tot_freqs * n_mice_w_dates * n_region_pairs);
end

for ix = 1:n_mice_w_dates
  next_data_dir = mice_all_info_table.Data_dir{ix};
  mouse_label = ['Mouse', mice_all_info_table.Animal_ID{ix}];
  date_label = mice_all_info_table.Date{ix};

  disp([ ...
    'Computing means for mouse ', ...
    mice_all_info_table.Animal_ID{ix}, ...
    ' on ', ...
    mice_all_info_table.Date{ix}, ...
    ' ...']);

  mouse_row_ax = ((ix - 1) * n_tot_freqs * n_region_pairs) + 1;
  mouse_row_bx = ix * n_tot_freqs * n_region_pairs;
  mouse_ixs = logical(zeros(size(mice_means_ts, 1), 1) );
  mouse_ixs(mouse_row_ax:mouse_row_bx) = true;

  mice_means_ts.mouse(mouse_ixs) = repelem({mouse_label}, n_tot_freqs * n_region_pairs, 1);
  mice_means_ts.date(mouse_ixs) = repelem({date_label}, n_tot_freqs * n_region_pairs, 1);

  mouse_ix_fname_ixs = ( ...
    strcmp(mice_filenames_table.Animal_ID, mice_all_info_table.Animal_ID{ix}) & ...
    strcmp(mice_filenames_table.Date, date_label) );

  mouse_ix_fnames_table = mice_filenames_table( ...
    mouse_ix_fname_ixs, ...
    [{'region_pair'}, {'filename'}]);

  mouse_ix_exp_phases_data = struct();

  for px = 1:n_exp_phases
    mouse_ix_exp_phases_data.(experiment_phases{px}) = NaN( ...
      n_phase_pts_vec(px), ...
      n_tot_freqs * n_region_pairs);
  end

  for jx = 1:n_region_pairs
    next_region_pair = all_region_pairs{jx};

    disp(['Computing means for region pair ', next_region_pair]);

    next_reg_pair_ixs = strcmp(mouse_ix_fnames_table.region_pair, next_region_pair);
    next_file_set = mouse_ix_fnames_table.filename(next_reg_pair_ixs);
    n_next_files = length(next_file_set);

    mouse_region_ixs = strcmp(mice_means_ts.region_pair(mouse_ixs), next_region_pair);

    freq_correctness_check = mice_means_ts.Hz(mouse_ixs);
    freq_correctness_check = freq_correctness_check(mouse_region_ixs);
    assert(all(freq_correctness_check == all_freqs') );

    selected_data = struct();

    for px = 1:n_exp_phases
      selected_data.(experiment_phases{px}) = NaN( ...
        n_phase_pts_vec(px), ...
        n_next_files * n_tot_freqs);
    end

    for fx = 1:n_next_files
      next_file = [next_data_dir, next_file_set{fx}];
      loaded_data = load(next_file);
      original_matrix = [];

      % Check if the specific variable name "CXY" exists.
      if isfield(loaded_data, 'CXY')
        original_matrix = loaded_data.CXY;
      else
        error('Matrix variable not found in the file: %s', next_file);
      end

      % Order the the data such that each group of n_next_files columns represent the coherence at
      % the same frequency in one Hz increments.
      col_ixs = fx:n_next_files:(n_tot_freqs * n_next_files);

      for px = 1:n_exp_phases
        exp_phase = experiment_phases{px};
        exp_phase_start = [exp_phase, '_start'];
        exp_phase_end = [exp_phase, '_end'];
        exp_phase_col_ax = mice_all_info_table.(exp_phase_start)(ix);
        exp_phase_col_bx = mice_all_info_table.(exp_phase_end)(ix);

        selected_data.(exp_phase)(:,col_ixs) = original_matrix( ...
          all_freqs, ...
          exp_phase_col_ax:exp_phase_col_bx)';
      end
    end

    for px = 1:n_exp_phases
      exp_phase = experiment_phases{px};
      mouse_ix_exp_phases_data.(exp_phase)(:,mouse_region_ixs) = reshape( ...
        mean(reshape(selected_data.(exp_phase)', n_next_files, []), 1, "omitnan"), ...
        n_tot_freqs, [])';
    end
  end

  for px = 1:n_exp_phases
    exp_phase = experiment_phases{px};
    exp_phases_data.(exp_phase)(:,mouse_ixs) = mouse_ix_exp_phases_data.(exp_phase);
    mouse_ix_exp_phases_data.(exp_phase) = [ ...
      mice_means_ts(mouse_ixs,:), ...
      array2table(mouse_ix_exp_phases_data.(exp_phase)')];
  end

  mouse_ix_save_filename = [save_dir, mouse_label, '_', date_label, '_mean_coherence_ts.mat'];
  save(mouse_ix_save_filename, '-struct', 'mouse_ix_exp_phases_data');
end

for px = 1:n_exp_phases
  exp_phase = experiment_phases{px};
  exp_phases_data.(exp_phase) = [ ...
    mice_means_ts, ...
    array2table(exp_phases_data.(exp_phase)')];
end

save_filename = [save_dir, 'mean_coherence_time_series_results.mat'];
save(save_filename, '-struct', 'exp_phases_data');

