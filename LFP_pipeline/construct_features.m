function construct_features( ...
  lfp_data_dir, ...
  save_dir, ...
  times_and_labels_table_fname, ...
  labels, ...
  max_freq, ...
  win_len, ...
  fs)

% Function call:
%   construct_features(lfp_data_dir, save_dir, times_and_labels_table_fname, labels, max_freq, win_len, fs);
%
% Description:
%   X
%
% Arguments:
%                   lfp_data_dir: Path to the averaged/ notch filtered LFP data files produced as a result
%                                 of calling the 'average_chans_and_filter_LFP' function.
%
%                       save_dir: Path to the directory where the features results files will be saved.
%
%   times_and_labels_table_fname: Filename of the CSV file containing the times and labels table. This table
%                                 should at least have columns: 'ID', 'date', 'start_time', 'end_time' and the
%                                 labels contained in the provided 'labels' variable. E.g.
%
%                                      ID   date Treatment Sex CGRP exp_phase start_time  end_time ...
%                                   MV211 240517   CGRP_A1   M    0    phase1          1       600 ...
%                                   MV211 240517   CGRP_A1   M    1    phase2        862      1019 ...
%                                     ...    ...       ... ...            ...        ...       ... ...
%                                   MV212 240517   CGRP_A1   M    0    phase1          1       600 ...
%                                   MV212 240517   CGRP_A1   M    1    phase2        847      1004 ...
%                                     ...    ...       ... ...            ...        ...       ... ...
%                                   MV214 240605   CGRP_A1   M    0    phase1          1       600 ...
%                                   MV214 240605   CGRP_A1   M    1    phase2        832       989 ...
%                                     ...    ...       ... ...            ...        ...       ... ...
%                                   MV217 240613   CGRP_A1   M    0    phase1          1       600 ...
%                                   MV217 240613   CGRP_A1   M    1    phase2        849      1006 ...
%                                     ...    ...       ... ...            ...        ...       ... ...
%                                   MV218 240621   CGRP_A1   M    0    phase1          1       600 ...
%                                   MV218 240621   CGRP_A1   M    1    phase2        815       972 ...
%                                     ...    ...       ... ...            ...        ...       ... ...
%                                   MV221 240621   CGRP_A1   F    0    phase1          1       600 ...
%                                   MV221 240621   CGRP_A1   F    1    phase2        831       988 ...
%                                     ...    ...       ... ...            ...        ...       ... ...
%
%                         labels: Vector of columns in the times and labels table used to label all time windows
%                                 within time segments indicated by the 'start_time' and 'end_time' columns.
%                                 For example, if the 'labels' variable is set to [{'CGRP'}, {'exp_phase'}],
%                                 and the above table is provided as the times and labels table, all time windows
%                                 for each mouse on each date would have corresponding labels determined by the
%                                 values in the "CGRP" and "exp_phase" columns.
%
%                       max_freq: Highest frequency (in Hz) of cross-power data to include.
%
%                        win_len: Integer indicating how many seconds of data to use in each window.
%
%                             fs: The sampling frequency of the LFP data.
%
% Example:
%   >>> lfp_data_dir = "./averaged_and_filtered_LFP_data/";
%   >>> save_dir = "./features_results/";
%   >>> times_and_labels_table_fname = "./times_and_labels_table.csv";
%   >>> labels = [{'CGRP'}, {'exp_phase'}];
%   >>> max_freq = 50;
%   >>> win_len = 1;
%   >>> fs = 1000;
%   >>> construct_features(lfp_data_dir, save_dir, times_and_labels_table_fname, labels, max_freq, win_len, fs);

  if isstring(max_freq) || ischar(max_freq)
    max_freq = str2num(max_freq);
  end

  if isstring(win_len) || ischar(win_len)
    win_len = str2num(win_len);
  end

  if isstring(fs) || ischar(fs)
    fs = str2num(fs);
  end

  % Check to make sure 'max_freq', 'win_len' and 'fs' are greater than zero.
  assert(max_freq > 0);
  assert(win_len > 0);
  assert(fs > 0);

  % Check to make sure 'win_len' and 'fs' are integers.
  assert(mod(win_len, 1) == 0);
  assert(mod(fs, 1) == 0);

  lfp_data_dir = char(lfp_data_dir);
  save_dir = char(save_dir);
  times_and_labels_table_fname = char(times_and_labels_table_fname);

  if ischar(labels)
    labels = [{labels}];
  end

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

  file_info = preprocess_LFP_support.extract_patterns_from_filename( ...
    lfp_filenames, ...
    re_patterns_table);

  % Make sure that the key variables for the LFP files are unique.
  lfp_group_ixs = findgroups(file_info(:,re_labels) );
  problem_lfp_group_ixs = find(groupcounts(lfp_group_ixs) > 1);

  if ~isempty(problem_lfp_group_ixs)
    problem_lfp_filenames = file_info.filename(ismember(lfp_group_ixs,  problem_lfp_group_ixs) );
    err_msg = [ ...
      'The following LFP files do not have unique key variables:\n', ...
      repmat('  %s\n', 1, numel(problem_lfp_filenames) )];

    error(sprintf(err_msg, problem_lfp_filenames{:}) );
  end

  table_opts = detectImportOptions(times_and_labels_table_fname);

  expected_colnames = [re_labels, labels, {'start_time'}, {'end_time'}];
  [~, ~, expected_ixs] = intersect(table_opts.VariableNames, expected_colnames);

  % Make sure that the provided table contains at least the variable names provided in the
  % 're_labels' variable defined above and the 'labels' variable passed to this function.
  if length(expected_ixs) ~= length(expected_colnames)
    problem_cols = expected_colnames(setdiff(1:length(expected_colnames), expected_ixs) );
    err_msg = [ ...
      'The following labels are missing from the provided times and labels table:\n', ...
      repmat('  %s\n', 1, numel(problem_cols) )];

    error(sprintf(err_msg, problem_cols{:}) );
  end

  ix_date = find(strcmp(table_opts.VariableNames, 'date') );

  % Check to make sure there aren't multiple date columns in the provided table.
  assert(length(ix_date) == 1);

  table_opts.VariableTypes{ix_date} = 'char';
  times_and_labels_table = readtable(times_and_labels_table_fname, table_opts);

  n_fft = win_len * fs;
  freq = [];
  n_freq = max_freq * win_len + 1;

  n_mice_groups = size(file_info, 1);

  for mx = 1:n_mice_groups
    mouse_group_mx = strjoin(table2cell(file_info(mx,re_labels) ), '_');

    disp(['Constructing features for group: ', mouse_group_mx, ' ...']);

    % Multiple checks are performed on the times provided by this table. First, all provided times
    % must be positive integers. Second, all end times must be at least as large as the corresponding
    % start times. Third, no time segments can overlap. If any of these three conditions are violated,
    % this group is skipped. Later, after the LFP data is loaded, if the times provided in the table
    % are larger than the number of seconds of data in the LFP data file, as many windows as possible
    % will be used. This means that some experimental phases could be left out, but the feature
    % construction will still proceed.
    times_and_labels_mx = innerjoin( ...
      file_info(mx,re_labels), ...
      times_and_labels_table, ...
      'Keys', ...
      re_labels);

    n_exp_phases_mx = size(times_and_labels_mx, 1);

    if n_exp_phases_mx == 0
      warn_msg = ['No times and labels table records for ', mouse_group_mx];
      warning(warn_msg);

      continue;
    end

    ixs_valid = ( ...
      (times_and_labels_mx.start_time > 0) & ...
      (mod(times_and_labels_mx.start_time, 1) == 0) & ...
      (times_and_labels_mx.end_time > 0) & ...
      (mod(times_and_labels_mx.end_time, 1) == 0) );

    if ~all(ixs_valid)
      warn_msg = [ ...
        'Invalid time segment(s) detected for %s.\n', ...
        'All start/end times must be positive integers.'];

      warning(sprintf(warn_msg, mouse_group_mx) );

      continue;
    end

    times_and_labels_mx = sortrows(times_and_labels_mx, 'start_time', 'ascend');

    ixs_valid = times_and_labels_mx.end_time >= times_and_labels_mx.start_time;

    if ~all(ixs_valid)
      ixs_problem_seg = find(~ixs_valid);
      warn_msg = [ ...
        'Invalid time segment(s) detected for %s.\nThe following ', ...
        'time segment(s) have end time less than corresponding start time:', ...
        repmat(' %d,', 1, length(ixs_problem_seg) )];

      warn_msg = warn_msg(1:(end-1) );

      warning( ...
        sprintf( ...
          warn_msg, ...
          mouse_group_mx, ...
          ixs_problem_seg(:) ));

      continue;
    end

    if n_exp_phases_mx > 1
      ixs_valid = times_and_labels_mx.end_time(1:(end-1) ) < times_and_labels_mx.start_time(2:end);

      if ~all(ixs_valid)
        ixs_problem_seg = find(~ixs_valid);
        warn_msg = [ ...
          'Invalid time segment(s) detected for %s.\nThe following ', ...
          'time segment(s) overlap(s) with another time segment:', ...
          repmat(' %d,', 1, length(ixs_problem_seg) )];

        warn_msg = warn_msg(1:(end-1) );

        warning( ...
          sprintf( ...
            warn_msg, ...
            mouse_group_mx, ...
            ixs_problem_seg(:) ));

        continue;
      end
    end

    filename_mx = [lfp_data_dir, file_info.filename{mx}];
    lfp_struct_mx = load(filename_mx);
    lfp_regions_mx = reshape(fieldnames(lfp_struct_mx), [], 1);
    n_regions_mx = numel(lfp_regions_mx);

    [regions_check, ixs_order] = sort(lfp_regions_mx);
    lfp_regions_mx = lfp_regions_mx(ixs_order);

    assert(all(strcmp(regions_check, lfp_regions_mx) ));

    lfp_mat_mx = cell2mat(struct2cell(lfp_struct_mx) );
    lfp_mat_mx = lfp_mat_mx(ixs_order,:);

    region_pairs_mx = nchoosek(lfp_regions_mx, 2);
    n_region_pairs_mx = size(region_pairs_mx, 1);

    region_pairs_mx = reshape( ...
      arrayfun( ...
        @(ix_row) strjoin(region_pairs_mx(ix_row,:), '_x_'), ...
        1:n_region_pairs_mx, ...
        'UniformOutput', false), ...
      [], 1);

    reg_pair_groups_ixs_table = [ ...
      cumsum([1, (n_regions_mx-1):-1:2]); ...
      cumsum((n_regions_mx-1):-1:1)]';

    n_time_max = floor(size(lfp_mat_mx, 2) / fs);
    ixs_valid = times_and_labels_mx.start_time <= n_time_max;

    if ~all(ixs_valid)
      ixs_problem_seg = find(~ixs_valid);
      warn_msg = [ ...
        'Invalid time segment(s) detected for %s.\nThere are ' ...
        '%d seconds of data in the LFP file. ', ...
        'The following time segments are listed as\n', ...
        'starting after this time and will be left out:', ...
        repmat(' %d,', 1, length(ixs_problem_seg) )];

      warn_msg = warn_msg(1:(end-1) );

      warning( ...
        sprintf( ...
          warn_msg, ...
          mouse_group_mx, ...
          n_time_max, ...
          ixs_problem_seg(:) ));
    end

    times_and_labels_mx = times_and_labels_mx(ixs_valid,:);
    n_exp_phases_mx = size(times_and_labels_mx, 1);

    if n_exp_phases_mx == 0
      warn_msg = ['No times and labels table records for ', mouse_group_mx];
      warning(warn_msg);

      continue;
    end

    ixs_valid = times_and_labels_mx.end_time <= n_time_max;

    if ~all(ixs_valid)
      ix_problem_seg = find(~ixs_valid);

      assert(ix_problem_seg == n_exp_phases_mx);

      warn_msg = [ ...
        'Invalid time segment detected for %s.\nThere are ' ...
        '%d seconds of data in the LFP file. The following ', ...
        'time segment is listed as\nending after this time ', ...
        'and will be shortened to end at %d seconds: %d'];

      warning( ...
        sprintf( ...
          warn_msg, ...
          mouse_group_mx, ...
          n_time_max, ...
          n_time_max, ...
          ix_problem_seg) );

      times_and_labels_mx.end_time(ix_problem_seg) = n_time_max;
    end

    n_sec_per_seg_mx = times_and_labels_mx.end_time - times_and_labels_mx.start_time + 1;
    n_win_per_seg_mx = reshape(floor(n_sec_per_seg_mx ./ win_len), [], 1);
    phase_win_ixs_table = [ ...
      cumsum([1; n_win_per_seg_mx(1:(end-1) )]), ...
      cumsum(n_win_per_seg_mx)];

    phase_win_ixs_table(n_win_per_seg_mx < 1,:) = nan;

    times_and_labels_mx.start_window = phase_win_ixs_table(:,1);
    times_and_labels_mx.end_window = phase_win_ixs_table(:,2);

    n_win_mx = sum(n_win_per_seg_mx);

    power_mx = nan([n_freq, n_regions_mx, n_win_mx]);
    fft_mx = nan([n_freq, n_regions_mx, n_win_mx]);
    coherence_mx = nan([n_freq, n_region_pairs_mx, n_win_mx]);
    labels_mx = array2table(cell([n_win_mx, numel(labels)]), 'VariableNames', labels);

    for px = 1:n_exp_phases_mx
      ix_time = times_and_labels_mx.start_time(px);

      ax_seg = phase_win_ixs_table(px,1);
      bx_seg = phase_win_ixs_table(px,2);
      n_win_mx_px = bx_seg - ax_seg + 1;

      if isnan(n_win_mx_px) || (n_win_mx_px < 1)
        continue;
      end

      update_msg = [ ...
        '  Constructing features for time windows %d through %d, ', ...
        'corresponding to time %d through %d seconds ...'];

      disp( ...
        sprintf( ...
          update_msg, ...
          ax_seg, ...
          bx_seg, ...
          times_and_labels_mx.start_time(px), ...
          times_and_labels_mx.end_time(px) ));

      for lx = 1:numel(labels)
        labels_mx(ax_seg:bx_seg,labels{lx}) = repelem(times_and_labels_mx.(labels{lx})(px), n_win_mx_px, 1);
      end

      % At each window, compute the fft, power and coherence for each brain region and each pair
      % of brain regions.
      for wx = 1:n_win_mx_px
        ix_win = ax_seg + wx - 1;

        ax_data = ((ix_time - 1) + (wx - 1) * win_len) * fs + 1;
        bx_data = ((ix_time - 1) + wx * win_len) * fs;

        % Slice LFP data to the next window.
        window_mat = lfp_mat_mx(:,ax_data:bx_data)';

        % Compute power.
        [power_mx_wx, freq_check] = cpsd(window_mat, window_mat, [], [], n_fft, fs);

        if length(freq) == 0
          freq = freq_check(freq_check <= max_freq);

          assert(length(freq) == n_freq);
        end

        assert(all(freq_check(1:n_freq) == freq) );

        %power_mx(:,:,ix_win) = power_mx_wx(1:n_freq,:) * 2 * pi / sqrt(fs);

        power_mx(:,:,ix_win) = power_mx_wx(1:n_freq,:);

        % Compute FFT.
        fft_mx_wx = fft(window_mat);
        fft_mx(:,:,ix_win) = fft_mx_wx(1:n_freq,:);

        % Compute coherence
        for rx = 1:(n_regions_mx - 1)
          ax_reg_pair = reg_pair_groups_ixs_table(rx,1);
          bx_reg_pair = reg_pair_groups_ixs_table(rx,2);

          assert((bx_reg_pair - ax_reg_pair + 1) == (n_regions_mx - rx) );

          [coherence_mx_wx, freq_check] = mscohere(window_mat(:,rx), window_mat(:,(rx+1):end), [], [], n_fft, fs);

          assert(all(freq_check(1:n_freq) == freq) );

          coherence_mx(:,ax_reg_pair:bx_reg_pair,ix_win) = coherence_mx_wx(1:n_freq,:);
        end
      end
    end

    save_struct = struct();
    save_struct.power = power_mx;
    save_struct.fft = fft_mx;
    save_struct.regions = lfp_regions_mx;
    save_struct.coherence = coherence_mx;
    save_struct.region_pairs = region_pairs_mx;
    save_struct.freq = freq;
    save_struct.n_win = n_win_mx;
    save_struct.fs = fs;
    save_struct.win_len = win_len;
    save_struct.group_ID = mouse_group_mx;

    for lx = 1:length(labels)
      save_struct.(labels{lx}) = labels_mx.(labels{lx});
    end

    save_struct.times_and_labels_table = times_and_labels_mx;

    save_filename = [save_dir, 'Mouse', mouse_group_mx, '_features.mat'];

    save(save_filename, '-struct', 'save_struct');
  end
end
