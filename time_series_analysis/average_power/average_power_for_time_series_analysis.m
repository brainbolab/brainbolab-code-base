% For each brain region provided, this script will compute the mean log-power time series over the
% indicated frequency ranges for each mouse on each date provided. The results are saved as a MATLAB
% struct with fields identifying tables containing mean log-power time series data in one Hz increments
% for different phases of an experiment.

% Change this to the path of the directory containing the mice data directories.
data_dir = [ ...
  '/Users', ...
  '/ikhultman', ...
  '/average_spectograms', ...
  '/test_data'];

% Change this to the path of the directory containing the mice CHANS files.
chans_dir = [ ...
  '/Users', ...
  '/ikhultman', ...
  '/average_spectograms', ...
  '/test_chans'];

% Change this to the path where the results should be saved.
save_dir = [ ...
  '/Users', ...
  '/ikhultman', ...
  '/average_spectograms', ...
  '/mean_time_series_results'];

% Change this to the path of the CSV file containing the mice IDs, dates and times at
% which the data will be trimmed.
times_table_filename = [ ...
  '/Users', ...
  '/ikhultman', ...
  '/average_spectograms', ...
  '/CGRPaba_Timepoints_230824.csv'];

% Change the values of this array such that the first element of each row is a mouse ID,
% and the second element of each row is the date corresponding to the data for that mouse
% for which the mean log-power will be computed. The mouse IDs in this array must be a subset
% of the mouse IDs provided in the times table.
mice_w_dates = [ ...
  {'M15'}, {'220204'}; ...
  {'M19'}, {'220811'}];

% Change this array to include the subset of brain regions for which the mean log-power will
% be computed. Make sure that the brain region labels included as part of the data filenames
% have the same format as those listed here.
brain_regions = { ...
  'BLA', ...
  'CeA', ...
  'LPBN'};

% Change this to include each frequency range in Hz over which to compute the mean log-power
% time series for each brain region for each mouse.
freq_ranges = [ ...
  {4:7}, ...
  {11:18}];

% Change this to include each phase of the experiment for which the mean log-power time series
% will be computed. Each phase label that appears here must have corresponding 'start' and 'end'
% labels in the provided times table. E.g. if one of the phase labels is 'Baseline', then there
% must be 'Baseline_start' and 'Baseline_end' fields in the times table identifying the first and
% last indices into each mouse's total dataset that correspond to that phase of the experiment.
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

if chans_dir(end) ~= file_sep
  chans_dir = [chans_dir, file_sep];
end

if save_dir(end) ~= file_sep
  save_dir = [save_dir, file_sep];
end

addpath(['.', file_sep, 'average_power_support']);

chans_filenames = ls_filenames_w_pattern(chans_dir, '.*CHANS.mat$');
chans_re_labels = [{'Animal_ID'}, {'Date'}];
chans_re_patterns_table = cell2table( ...
  [{'(?i)(?<=^mouse)[^_]*'}, {'(?<=_)[0-9]{6}(?=_)'}], ...
  'VariableNames', chans_re_labels);

chans_labels_info = extract_patterns_from_filename( ...
  chans_filenames, ...
  chans_re_patterns_table);

chans_labels_info = renamevars(chans_labels_info, 'filename', 'chans_filename');

exp_phase_colnames = cellfun( ...
  @(exp_phase) [{[exp_phase, '_start']}, {[exp_phase, '_end']}], ...
  experiment_phases, ...
  'UniformOutput', false);

exp_phase_colnames = [exp_phase_colnames{:}];
expected_colnames = [chans_re_labels, exp_phase_colnames];

table_opts = detectImportOptions(times_table_filename);
[~, ~, expected_ixs] = intersect(table_opts.VariableNames, expected_colnames);

if isrow(expected_ixs)
  expected_ixs = expected_ixs';
end

% Make sure that the provided table contains the variable names provided in the
% expected_colnames variable defined above. It's fine if the table contains more
% than just the variable names provided (and order doesn't matter), but it must
% contain at least those expected variable names.
assert( ...
  length(expected_ixs) == length(expected_colnames) && ...
  all(sort(expected_ixs) == (1:length(expected_colnames) )') );

date_ix = find(strcmp(table_opts.VariableNames, 'Date') );

% Check to make sure there aren't multiple date columns in the provided table.
assert(length(date_ix) == 1);

table_opts.VariableTypes{date_ix} = 'char';
times_table = readtable(times_table_filename, table_opts);

mice_w_dates_table = join( ...
  cell2table(mice_w_dates, 'VariableNames', chans_re_labels), ...
  chans_labels_info);

n_mice_w_dates = size(mice_w_dates_table, 1);
mice_w_dates_table.Data_dir = cell(n_mice_w_dates, 1);

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
    'Spectogram', ...
    file_sep];

  mouse_ix_fnames = ls_filenames_w_pattern( ...
    mice_w_dates_table.Data_dir{ix}, ...
    [regions_re_pattern, '.*\.mat$']);

  mouse_ix_chans_re_table = cell2table( ...
    [{[regions_re_pattern, '_[0-9]*']}, {regions_re_pattern}], ...
    'VariableNames', [{'channel'}, {'region'}]);

  mouse_ix_fnames_w_chans_table = extract_patterns_from_filename( ...
    mouse_ix_fnames, ...
    mouse_ix_chans_re_table);

  mouse_ix_chans_info = load([chans_dir, mice_w_dates_table.chans_filename{ix}]);

  if isrow(mouse_ix_chans_info.CHANNAMES)
    mouse_ix_chans_info.CHANNAMES = mouse_ix_chans_info.CHANNAMES';
  end

  if isrow(mouse_ix_chans_info.CHANACTIVE)
    mouse_ix_chans_info.CHANACTIVE = mouse_ix_chans_info.CHANACTIVE';
  end

  mouse_ix_chans_ixs = ~cellfun( ...
    @isempty, ...
    regexp(mouse_ix_chans_info.CHANNAMES, regions_re_pattern) );

  if isrow(mouse_ix_chans_ixs)
    mouse_ix_chans_ixs = mouse_ix_chans_ixs';
  end

  mouse_ix_chans_ixs = mouse_ix_chans_ixs & logical(mouse_ix_chans_info.CHANACTIVE);
  mouse_ix_chans_table = table( ...
    mouse_ix_chans_info.CHANNAMES(mouse_ix_chans_ixs), ...
    'VariableNames', {'channel'});

  mouse_ix_fnames_table = join(mouse_ix_chans_table, mouse_ix_fnames_w_chans_table);
  n_mouse_ix_fnames = size(mouse_ix_fnames_table, 1);

  mouse_ix_fnames_table.Animal_ID = repelem( ...
    mice_w_dates_table.Animal_ID(ix), n_mouse_ix_fnames, 1);

  mouse_ix_fnames_table.Date = repelem( ...
    mice_w_dates_table.Date(ix), n_mouse_ix_fnames, 1);

  mouse_ix_fnames_table_cols = [{'Animal_ID'}, {'Date'}, {'region'}, {'filename'}];
  mouse_ix_fnames_table = mouse_ix_fnames_table(:,mouse_ix_fnames_table_cols);

  if isempty(mice_filenames_table)
    mice_filenames_table = mouse_ix_fnames_table;
  else
    mice_filenames_table = [mice_filenames_table; mouse_ix_fnames_table];
  end
end

mice_all_info_table = join(mice_w_dates_table, times_table);

n_regions = length(brain_regions);
n_freq_ranges = length(freq_ranges);
n_freqs_each = arrayfun(@(rx) length(freq_ranges{rx}), 1:n_freq_ranges);

all_freqs = [freq_ranges{:}];
n_tot_freqs = length(all_freqs);

mice_means_ts_colnames = [ {'mouse'},   {'date'}, {'freq_band'},    {'Hz'}, {'region'}];
mice_means_ts_vartypes = [{'string'}, {'string'},     {'int64'}, {'int64'}, {'string'}];

assert(numel(mice_means_ts_colnames) == numel(mice_means_ts_vartypes) );

mice_means_ts = table( ...
  'Size', [n_tot_freqs * n_mice_w_dates * n_regions, numel(mice_means_ts_colnames)], ...
  'VariableTypes', mice_means_ts_vartypes, ...
  'VariableNames', mice_means_ts_colnames);

frq_grp_ids = arrayfun( ...
  @(rx) repelem(rx, n_freqs_each(rx) ), ...
  1:n_freq_ranges, ...
  'UniformOutput', false);

frq_grp_ids = repmat(repmat([frq_grp_ids{:}], 1, n_regions), 1, n_mice_w_dates);
mice_means_ts.freq_band = frq_grp_ids';
mice_means_ts.Hz = repmat(all_freqs', n_regions * n_mice_w_dates, 1);
mice_means_ts.region = repmat( ...
  repelem(brain_regions, 1, n_tot_freqs)', n_mice_w_dates, 1);

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
    n_tot_freqs * n_mice_w_dates * n_regions);
end

for ix = 1:n_mice_w_dates
  next_data_dir = mice_all_info_table.Data_dir{ix};
  mouse_label = ['Mouse', mice_all_info_table.Animal_ID{ix}];
  date_label = mice_all_info_table.Date{ix};

  disp([ ...
    'Computing log-power time series means for mouse ', ...
    mice_all_info_table.Animal_ID{ix}, ...
    ' on ', ...
    mice_all_info_table.Date{ix}, ...
    ' ...']);

  mouse_row_ax = ((ix - 1) * n_tot_freqs * n_regions) + 1;
  mouse_row_bx = ix * n_tot_freqs * n_regions;
  mouse_ixs = logical(zeros(size(mice_means_ts, 1), 1) );
  mouse_ixs(mouse_row_ax:mouse_row_bx) = true;

  mice_means_ts.mouse(mouse_ixs) = repelem({mouse_label}, n_tot_freqs * n_regions, 1);
  mice_means_ts.date(mouse_ixs) = repelem({date_label}, n_tot_freqs * n_regions, 1);

  mouse_ix_fname_ixs = ( ...
    strcmp(mice_filenames_table.Animal_ID, mice_all_info_table.Animal_ID{ix}) & ...
    strcmp(mice_filenames_table.Date, date_label) );

  mouse_ix_fnames_table = mice_filenames_table( ...
    mouse_ix_fname_ixs, ...
    [{'region'}, {'filename'}]);

  mouse_ix_exp_phases_data = struct();

  for px = 1:n_exp_phases
    mouse_ix_exp_phases_data.(experiment_phases{px}) = NaN( ...
      n_phase_pts_vec(px), ...
      n_tot_freqs * n_regions);
  end

  for jx = 1:n_regions
    next_region = brain_regions{jx};

    disp(['Computing means for region ', next_region, ' ...']);

    next_reg_ixs = strcmp(mouse_ix_fnames_table.region, next_region);
    next_file_set = mouse_ix_fnames_table.filename(next_reg_ixs);
    n_next_files = length(next_file_set);

    mouse_region_ixs = strcmp(mice_means_ts.region(mouse_ixs), next_region);

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
      loaded_data = load(next_file, 'PowerVar');
      original_matrix = [];

      % Check if the specific variable name "PowerVar.PowerdB" exists.
      if isfield(loaded_data, 'PowerVar') && isfield(loaded_data.PowerVar, 'PowerdB')
        original_matrix = loaded_data.PowerVar.PowerdB;
      else
        error('Matrix variable not found in the file: %s', next_file);
      end

      % Order the the data such that each group of n_next_files columns represent the power at
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
      neg_pow_ixs = selected_data.(exp_phase) <= 0;
      selected_data.(exp_phase)(neg_pow_ixs) = NaN;
      selected_data.(exp_phase) = log(selected_data.(exp_phase) );

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

  mouse_ix_save_filename = [save_dir, mouse_label, '_', date_label, '_mean_logpower_ts.mat'];
  save(mouse_ix_save_filename, '-struct', 'mouse_ix_exp_phases_data');
end

for px = 1:n_exp_phases
  exp_phase = experiment_phases{px};
  exp_phases_data.(exp_phase) = [ ...
    mice_means_ts, ...
    array2table(exp_phases_data.(exp_phase)')];
end

save_filename = [save_dir, 'mean_logpower_time_series_results.mat'];
save(save_filename, '-struct', 'exp_phases_data');

