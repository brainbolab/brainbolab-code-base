function compute_coherence_for_averaged_LFP( ...
  lfp_data_dir, ...
  save_dir, ...
  max_plot_freq, ...
  win_len, ...
  win_step, ...
  fs)
% Function call:
%   compute_coherence_for_averaged_LFP(lfp_data_dir, save_dir, max_plot_freq, win_len, win_step, fs)
%
% Description:
%   Compute the magnitude-squared coherence between each pair of regions at each time window in each LFP file
%   located in the 'lfp_data_dir' directory. The coherence results and figures will be saved in the
%   subdirectories of the 'save_dir' directory. The 'win_len' and 'win_step' variables are used to indicate
%   window length (in seconds) and the amount of overlap (in seconds) between subsequent time windows. If
%   you do not want windows to overlap, set the 'win_len' and 'win_step' variables to be the same.
%
% Arguments:
%    lfp_data_dir: Path to the averaged/ notch filtered LFP data files produced as a result of calling
%                  the 'average_chans_and_filter_LFP' function.
%
%        save_dir: Path to the directory where the coherence results files will be saved.
%
%   max_plot_freq: Highest frequency (in Hz) of coherence data to include in figures.
%
%         win_len: Integer indicating how many seconds of data to use in each window. This should be at least
%                  as large as 'win_step'.
%
%        win_step: Integer indicating the step size (in seconds) between starting times of subsequent windows.
%                  This should be less than or equal to 'win_len'.
%
%              fs: The sampling frequency of the LFP data.
%
% Example:
%   >>> lfp_data_dir = "./averaged_and_filtered_LFP_data/";
%   >>> save_dir = "./coherence_results/";
%   >>> max_plot_freq = 50;
%   >>> win_len = 1;
%   >>> win_step = 1;
%   >>> fs = 1000;
%   >>> compute_coherence_for_averaged_LFP(lfp_data_dir, save_dir, max_plot_freq, win_len, win_step, fs);

  % Check to make sure 'max_plot_freq', 'win_len', 'win_step' and 'fs' are greater than zero.
  assert(max_plot_freq > 0);
  assert(win_len > 0);
  assert(win_step > 0);
  assert(fs > 0);

  % Check to make sure 'win_len', 'win_step' and 'fs' are integers.
  assert(mod(win_len, 1) == 0);
  assert(mod(win_step, 1) == 0);
  assert(mod(fs, 1) == 0);

  if (win_len < win_step)
    warning_msg = 'Variable win_len = %d which is less than variable win_step = %d. This will work but does not make sense.';
    warning(sprintf(warning_msg, win_len, win_step) );
  end

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

  % Create results directories.
  coherence_save_dir = [save_dir, 'coherence_results', file_sep];

  if ~isfolder(coherence_save_dir)
    mkdir(coherence_save_dir);

    assert(isfolder(coherence_save_dir) );
  end

  plots_save_dir = [save_dir, 'figures', file_sep];

  if ~isfolder(plots_save_dir)
    mkdir(plots_save_dir);

    assert(isfolder(plots_save_dir) );
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

  for fx = 1:size(labels_info, 1)
    disp_str = sprintf( ...
      'Computing coherence for Mouse%s on %s ...', ...
      labels_info.ID{fx}, ...
      labels_info.date{fx});

    disp(disp_str);

    lfp_filename_fx = [lfp_data_dir, labels_info.filename{fx}];
    lfp_struct = load(lfp_filename_fx);

    regions = fieldnames(lfp_struct);
    n_regions = numel(regions);
    n_region_pairs = nchoosek(n_regions, 2);

    n_data_pts = length(lfp_struct.(regions{1}) );
    n_blocks = floor(n_data_pts / fs);
    n_windows = floor((n_blocks - win_len) / win_step) + 1;
    nfft = win_len * fs;
    n_freq = floor(nfft / 2) + 1;

    freq = [];
    region_pairs = cell(n_region_pairs, 1);
    win_times = [];

    CXY = nan(n_freq, n_region_pairs, n_windows);
    ix_start = 1;

    % Compute coherence between each pair of regions at each time window.
    for ix_in = 1:(n_regions - 1)
      ixs_out = (ix_in + 1):n_regions;
      in_region = regions{ix_in};
      out_regions = regions(ixs_out);
      ix_end = ix_start + numel(out_regions) - 1;

      [CXY_ix, freq_ix, region_labels, win_times_ix] = preprocess_LFP_support.compute_coherence( ...
        lfp_struct, ...
        in_region, ...
        out_regions, ...
        win_len, ...
        win_step, ...
        fs);

      if isempty(freq)
        freq = reshape(freq_ix, [], 1);
      end

      if isempty(win_times)
        win_times = reshape(win_times_ix, [], 1);
      end

      CXY(:,ix_start:ix_end,:) = CXY_ix;

      if ischar(region_labels)
        region_labels = {region_labels};
      end

      region_pairs(ix_start:ix_end) = reshape(region_labels, [], 1);
      ix_start = ix_end + 1;
    end

    disp('Saving coherence results and plotting figures ...');

    save_filename = sprintf( ...
      '%sMouse%s_%s_coherence_results.mat', ...
      coherence_save_dir, ...
      labels_info.ID{fx}, ...
      labels_info.date{fx});

    % Save the coherence data.
    save(save_filename, 'CXY', 'freq', 'region_pairs', 'win_times');

    plots_save_dir_fx = sprintf( ...
      ['%sMouse%s', file_sep, '%s', file_sep], ...
      plots_save_dir, ...
      labels_info.ID{fx}, ...
      labels_info.date{fx});

    if ~isfolder(plots_save_dir_fx)
      mkdir(plots_save_dir_fx);

      assert(isfolder(plots_save_dir_fx) );
    end

    ixs_plot_freq = find(freq <= max_plot_freq);
    n_plot_freq = length(ixs_plot_freq);

    % Plot the coherence spectograms for each pair of regions.
    for rx = 1:n_region_pairs
      CXY_rx = reshape(CXY(ixs_plot_freq,rx,:), n_plot_freq, n_windows);
      fig_title = sprintf( ...
        'Mouse%s %s %s', ...
        labels_info.ID{fx}, ...
        labels_info.date{fx}, ...
        char(join(split(region_pairs{rx}, '_x_'), ':') ));

      plot_rx = figure();
      pcolor(win_times, freq(ixs_plot_freq), CXY_rx);
      shading flat;
      cb = colorbar();
      cb.Label.String = 'Coherence';
      title(fig_title, 'Interpreter', 'none');
      xlabel('Time (seconds)');
      ylabel('Frequency (Hz)');
      plot_save_filename = [plots_save_dir_fx, region_pairs{rx}, '.bmp'];
      saveas(plot_rx, plot_save_filename, 'bmp');
      close(plot_rx);
    end
  end
end
