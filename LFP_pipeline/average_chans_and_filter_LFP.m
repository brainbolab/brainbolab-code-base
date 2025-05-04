function average_chans_and_filter_LFP(lfp_data_dir, chans_dir, save_dir, fs)
% Function call:
%   average_chans_and_filter_LFP(lfp_data_dir, chans_dir, save_dir, fs)
%
% Description:
%   For each LFP file in the provided data directory, this function computes the average LFP signal
%   across the channels considered "good" by the corresponding CHANS file for each brain region and
%   then notch filters the averaged LFP signals at 60 Hz. The averaged and filtered LFP data files
%   are then saved to the provided directory.
%
% Arguments:
%   lfp_data_dir: Path to the original LFP data files.
%
%      chans_dir: Path to the corresponding CHANS files. Each LFP file must have a corresponding
%                 CHANS file or it will be ignored.
%
%       save_dir: Path to the directory where the averaged and notch filtered LFP data files will
%                 be saved.
%
%             fs: The sampling frequency of the LFP data.
%
% Example:
%   >>> lfp_data_dir = "./LFP_data/";
%   >>> chans_dir = "./CHANS/";
%   >>> save_dir = "./averaged_and_filtered_LFP_data/";
%   >>> fs = 1000;
%   >>> average_chans_and_filter_LFP(lfp_data_dir, chans_dir, save_dir, fs);

  lfp_data_dir = char(lfp_data_dir);
  chans_dir = char(chans_dir);
  save_dir = char(save_dir);

  file_sep = '/';

  if ispc
    file_sep = '\';
  end

  if (lfp_data_dir(end) ~= file_sep)
    lfp_data_dir = [lfp_data_dir, file_sep];
  end

  if (chans_dir(end) ~= file_sep)
    chans_dir = [chans_dir, file_sep];
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

  chans_filenames = preprocess_LFP_support.ls_filenames_w_pattern(chans_dir, '.*CHANS.*\.mat$');
  chans_filenames = sort(chans_filenames);

  lfp_filenames = reshape(lfp_filenames, numel(lfp_filenames), 1);
  chans_filenames = reshape(chans_filenames, numel(chans_filenames), 1);

  %% You may need to change this table depending on what information is contained in the filenames.
  re_labels = [{'ID'}, {'date'}];
  re_patterns_table = cell2table( ...
    [{'(?i)(?<=^mouse)[^_]*'}, {'(?<=_)[0-9]{6}(?=_)'}], ...
    'VariableNames', re_labels);

  labels_info = preprocess_LFP_support.extract_patterns_from_filename( ...
    lfp_filenames, ...
    re_patterns_table);

  chans_labels_info = preprocess_LFP_support.extract_patterns_from_filename( ...
    chans_filenames, ...
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

  % Make sure that the key variables for the CHANS files are unique.
  chans_group_ixs = findgroups(chans_labels_info(:,re_labels) );
  problem_chans_group_ixs = find(groupcounts(chans_group_ixs) > 1);

  if ~isempty(problem_chans_group_ixs)
    problem_chans_filenames = chans_labels_info.filename(ismember(chans_group_ixs,  problem_chans_group_ixs) );
    err_msg = [ ...
      'The following CHANS files do not have unique key variables:\n', ...
      repmat('  %s\n', 1, numel(problem_chans_filenames) )];

    error(sprintf(err_msg, problem_chans_filenames{:}) );
  end

  % Make sure that all of the ID's and dates for the input LFP files are the same ID's and
  % dates for the input CHANS files.
  [~, diff_ixs] = setdiff( ...
    labels_info(:,re_labels), ...
    chans_labels_info(:,re_labels) );

  if ~isempty(diff_ixs)
    problem_lfp_filenames = labels_info.filename(diff_ixs);
    warning_msg = [ ...
      'The following LFP files have no matching CHANS files and will be left out:\n', ...
      repmat('  %s\n', 1, length(diff_ixs) )];

    warning(sprintf(warning_msg, problem_lfp_filenames{:}) );
  end

  [~, diff_ixs] = setdiff( ...
    chans_labels_info(:,re_labels), ...
    labels_info(:,re_labels) );

  if ~isempty(diff_ixs)
    problem_chans_filenames = chans_labels_info.filename(diff_ixs);
    warning_msg = [ ...
      'The following CHANS files have no matching LFP files and will be left out:\n', ...
      repmat('  %s\n', 1, length(diff_ixs) )];

    warning(sprintf(warning_msg, problem_chans_filenames{:}) );
  end

  labels_info = renamevars(labels_info, 'filename', 'lfp_filename');
  chans_labels_info = renamevars(chans_labels_info, 'filename', 'chans_filename');

  labels_info = innerjoin(labels_info, chans_labels_info);

  n_files = size(labels_info, 1);

  expected_chans_fields = [{'CHANNAMES'}, {'CHANACTIVE'}];

  for fx = 1:n_files
    lfp_filename_fx = [lfp_data_dir, labels_info.lfp_filename{fx}];
    chans_filename_fx = [chans_dir, labels_info.chans_filename{fx}];

    lfp_struct = load(lfp_filename_fx);
    chans_info = load(chans_filename_fx);

    % Make sure that the 'chans_info' struct has the correct fields.
    assert(isempty(setdiff(fieldnames(chans_info), expected_chans_fields) ));
    assert(isempty(setdiff(expected_chans_fields, fieldnames(chans_info) )));

    % Make sure that the 'chans_info' and 'lfp_struct' structs represent the same channels.
    assert(isempty(setdiff(fieldnames(lfp_struct), chans_info.CHANNAMES) ));
    assert(isempty(setdiff(chans_info.CHANNAMES, fieldnames(lfp_struct) )));

    % Obtain the mean signals over all channels corresponding to each brain region.
    lfp_struct = preprocess_LFP_support.average_LFP(lfp_struct, chans_info);

    % Filter out 60Hz signal.
    lfp_struct = preprocess_LFP_support.notch_filter60_LFP(lfp_struct, fs);

    [~, save_filename_fx] = fileparts(labels_info.lfp_filename{fx});
    save_filename_fx = [save_dir, save_filename_fx, '_averaged_and_notch_filtered.mat'];
    save(save_filename_fx, '-struct', 'lfp_struct');
  end
end
