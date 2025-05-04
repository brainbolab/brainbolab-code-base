function compute_directionality_for_averaged_LFP( ...
  lfp_data_dir, ...
  save_dir, ...
  fs)
% Function call:
%   compute_directionality_for_averaged_LFP(lfp_data_dir, save_dir, fs)
%
% Description:
%   X
%
% Arguments:
%   lfp_data_dir: Path to the averaged/ notch filtered LFP data files produced as a result of calling
%                 the 'average_chans_and_filter_LFP' function.
%
%       save_dir: Path to the directory where the coherence results files will be saved.
%
%             fs: The sampling frequency of the LFP data.
% Example:
%   >>> lfp_data_dir = "./averaged_and_filtered_LFP_data/";
%   >>> save_dir = "./coherence_results/";
%   >>> fs = 1000;
%   >>> compute_directionality_for_averaged_LFP(lfp_data_dir, save_dir, fs);

  addpath(genpath('./directionality_support') );

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
    [{'(?i)(?<=^mouse)[^_]*'}, {'(?<=_)[0-9]{6}(?=_)'}], ...
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
      'Computing directionality for Mouse%s on %s ...', ...
      labels_info.ID{fx}, ...
      labels_info.date{fx});

    disp(disp_str);

    lfp_filename_fx = [lfp_data_dir, labels_info.filename{fx}];
    lfp_struct = load(lfp_filename_fx);

    regions = fieldnames(lfp_struct);
    n_regions = numel(regions);
    n_region_pairs = nchoosek(n_regions, 2);

  end

  % loc_1 and loc_2 provide the character indices into the file path 
  % string demarcating the actual file name within the full path.
  loc_1 = strfind(lfp_filename, sep);
  if isempty(loc_1)
    loc_1 = 0;
  end
  loc_1 = loc_1(end) + 1;

  loc_2 = regexp(lfp_filename, '\.');
  if isempty(loc_2)
    loc_2 = length(lfp_filename) + 1;
  end
  loc_2 = loc_2(1) - 1;

  lfp_struct = open(lfp_filename);

  disp(['LFP data loaded from ' lfp_filename]);

  fs_range = [1 50];
  time_ranges = [0 600]; % in seconds
  FS = 1000;

  if ~isempty(time_ranges)
    lfp_slices = {};
    for ix = 1:size(time_ranges, 1)
      disp(['Taking data from ' num2str(time_ranges(ix,1) ) 's to '...
            num2str(time_ranges(ix,2) ) 's ...']);

      ax = (time_ranges(ix,1) * FS) + 1;
      bx = (time_ranges(ix,2) * FS);
      lfp_slices{ix} = split_and_order_lfp(lfp_struct, ax, bx);
    end

    if (numel(lfp_slices) == 1)
      lfp_struct = lfp_slices{1};
    else
      lfp_struct = combine_lfp_structs(lfp_slices);
    end
  end

  brain_regions = regexprep(fieldnames(lfp_struct), "_[0-9]+$", "");
  brain_regions = unique(brain_regions(:,1) );

  assert(numel(brain_regions) > 1, ...
         ['ERROR in run_phase_offset()' newline ...
          'Only ' numel(brain_regions) ' brain regions found in the ' ...
          'provided LFP data file; at least two brain regions required.']);

  % Compare all of the brain regions.
  for rx_1 = 1:length(brain_regions)
    region_1 = brain_regions{rx_1};

    for rx_2 = (rx_1+1):length(brain_regions)
      region_2 = brain_regions{rx_2};

      disp(['Getting phase offsets for regions ' region_1 ...
            ' and ' region_2 ' ...']);

      % Get phase offsets for every pair of channels for the given regions
      % and time the process.
      tic;
      phase_offsets = get_phase_offsets( ...
        lfp_struct, ...
        region_1, ...
        region_2, ...
        fs_range);
      toc;

      save_filename = [save_dir lfp_filename(loc_1:loc_2) ...
                        '_' region_1 '_' region_2 '.mat'];

      disp(['Saving results to ' save_filename]);

      save(save_filename, 'phase_offsets');
    end
  end
end
