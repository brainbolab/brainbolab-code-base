function trim_LFP_from_tables( ...
  lfp_data_dir, ...
  save_dir, ...
  times_table_filename, ...
  fs)
% Function call:
%   trim_LFP_from_tables(lfp_data_dir, save_dir, times_table_filename)
%
% Description:
%   X
%
% Arguments:
%           lfp_data_dir: Path to the averaged/ notch filtered LFP data files produced as a result
%                         of calling the 'average_chans_and_filter_LFP' function.
%
%               save_dir: Path to the directory where the coherence results files will be saved.
%
%   times_table_filename:
%
%                     fs: The sampling frequency of the LFP data.
%
% Example:
%   >>> lfp_data_dir = "./averaged_and_filtered_LFP_data/";
%   >>> save_dir = "./trimmed_LFP_results/";
%   >>> times_table_filename = "./times_table.csv";
%   >>> trim_LFP_from_tables(lfp_data_dir, save_dir, times_table_filename);

  % Check to make sure 'fs' is greater than zero.
  assert(fs > 0);

  % Check to make sure 'fs' is an integer.
  assert(mod(fs, 1) == 0);

  lfp_data_dir = char(lfp_data_dir);
  save_dir = char(save_dir);

  file_sep = '/';

  if ispc
    file_sep = '\';
  end

  if (lfp_data_dir(end) ~= file_sep)
    lfp_data_dir = [lfp_data_dir, file_sep];
  end

  if (save_dir(end) ~= file_sep)
    save_dir = [save_dir, file_sep];
  end

  if ~isfolder(save_dir)
    mkdir(save_dir);

    assert(isfolder(save_dir) );
  end

  lfp_filenames = preprocess_LFP_support.ls_filenames_w_pattern(lfp_data_dir, '.*LFP.*\.mat$');
  lfp_filenames = sort(lfp_filenames);

  %% You may need to change this table depending on what information is contained in the filenames.
  re_labels = [{'ID'}, {'date'}];
  re_patterns_table = cell2table( ...
    [{'(?i)(?<=(^mouse)|(_mouse))[^_]*'}, {'(?<=_|^)[0-9]{6}(?=_)'}], ...
    'VariableNames', re_labels);

  labels_info = preprocess_LFP_support.extract_patterns_from_filename( ...
    lfp_filenames, ...
    re_patterns_table);

  % Make sure that the key variables for the LFP files are unique.
  lfp_group_ixs = findgroups(labels_info(:,re_labels) );
  problem_lfp_group_ixs = find(groupcounts(lfp_group_ixs) > 1);

  if ~isempty(problem_lfp_group_ixs)
    problem_lfp_filenames = labels_info.filename(ismember(lfp_group_ixs,  problem_lfp_group_ixs) );
    err_msg = [ ...
      'The following LFP files do not have unique key variables:\n', ...
      repmat('  %s\n', 1, numel(problem_lfp_filenames) )];

    error(sprintf(err_msg, problem_lfp_filenames{:}) );
  end

  key_colnames = [{'ID'}, {'date'}];
  table_opts = detectImportOptions(times_table_filename);

  [~, ~, expected_ixs] = intersect(table_opts.VariableNames, key_colnames);

  % Make sure that the provided table contains at least the variable names provided in the
  % key_colnames variable defined above.
  assert(length(expected_ixs) == length(key_colnames) );

  ix_date = find(strcmp(table_opts.VariableNames, 'date') );

  % Check to make sure there aren't multiple date columns in the provided table.
  assert(length(ix_date) == 1);

  table_opts.VariableTypes{ix_date} = 'char';
  times_table = readtable(times_table_filename, table_opts);

  exp_phases_present = cellfun( ...
    @(var_ix) char(regexp(var_ix, '.*(?=_((start)|(end))$)', 'match') ), ...
    times_table.Properties.VariableNames, ...
    'UniformOutput', false);

  ixs_exp_phases = cellfun( ...
    @(var_ix) length(var_ix) > 0, ...
    exp_phases_present);

  [grp_counts, exp_phases] = groupcounts(reshape(exp_phases_present(ixs_exp_phases), [], 1) );
  exp_phases = exp_phases(grp_counts == 2);
  ixs_exp_phases = false(length(exp_phases), 1);

  for px = 1:length(exp_phases)
    ixs_exp_phases(px) = ( ...
      (sum(strcmp(times_table.Properties.VariableNames, [exp_phases{px}, '_start']) ) == 1) & ...
      (sum(strcmp(times_table.Properties.VariableNames, [exp_phases{px}, '_end']) ) == 1) );
  end

  exp_phases = exp_phases(ixs_exp_phases);
  n_exp_phases = numel(exp_phases);

  assert(n_exp_phases > 0);

  for px = 1:n_exp_phases
    save_dir_px = [save_dir, exp_phases{px}];

    if ~isfolder(save_dir_px)
      mkdir(save_dir_px);

      assert(isfolder(save_dir_px) );
    end
  end

  times_w_labels_info = join(labels_info, times_table, 'Keys', key_colnames);
  n_filenames = size(times_w_labels_info, 1);

  for fx = 1:n_filenames
    filename_fx = [lfp_data_dir, times_w_labels_info.filename{fx}];

    disp(['Trimming data from file: ', filename_fx]);

    lfp_struct = load(filename_fx);
    lfp_regions = fieldnames(lfp_struct);
    n_regions = numel(lfp_regions);
    lfp_mat = cell2mat(struct2cell(lfp_struct) );
    n_data_pts = size(lfp_mat, 2);

    [~, save_filename_fx] = fileparts(times_w_labels_info.filename{fx});

    correctness_check = false(1, n_exp_phases);

    for px = 1:n_exp_phases
      save_filename_fx_px = [ ...
        save_dir, ...
        exp_phases{px}, ...
        file_sep, ...
        save_filename_fx, ...
        '_', ...
        exp_phases{px}, ...
        '.mat'];

      ix_start = (times_w_labels_info.([exp_phases{px}, '_start'])(fx) - 1) * fs + 1;
      ix_end = times_w_labels_info.([exp_phases{px}, '_end'])(fx) * fs;

      if (ix_end <= ix_start)
        warning_msg = [ ...
          'The end time must be greater than the start time for phase ', ...
          exp_phases{px}, ...
          ' in file: ', ...
          filename_fx];

        warning(warning_msg);

        continue
      end

      if (ix_end > n_data_pts)
        warning_msg = sprintf( ...
          ['Trying to index beyond the total number of data points in file %s\n' ...
          'There are %f seconds of data. Trying to index up to %f seconds.'], ...
          filename_fx, ...
          n_data_pts / fs, ...
          times_w_labels_info.([exp_phases{px}, '_end'])(fx) );

        warning(warning_msg);

        continue;
      end

      lfp_struct_px = cell2struct( ...
        mat2cell( ...
          lfp_mat(:,ix_start:ix_end), ...
          ones(1, n_regions) ), ...
        lfp_regions);

      correctness_check_px = false(1, n_regions);

      for rx = 1:n_regions
        reg_rx = lfp_regions{rx};
        correctness_check_px(rx) = all(lfp_struct.(reg_rx)(ix_start:ix_end) == lfp_struct_px.(reg_rx) );
      end

      correctness_check(px) = all(correctness_check_px);

      save(save_filename_fx_px, '-struct', 'lfp_struct_px');
    end

    disp(all(correctness_check) );
  end
end
