function fname_info = extract_patterns_from_filename( ...
  filenames, ...
  re_patterns_table)

  if isa(filenames, 'char')
    filenames = {filenames};
  end

  filenames = reshape(filenames, numel(filenames), 1);

  sep = '/';
  if ispc
    sep = '\';
  end

  % Find the indices in the provided filenames that separate the paths 
  % from the filenames.
  fname_ixs = strfind(filenames, sep);
  fname_ixs(cellfun(@(ixs) isempty(ixs), fname_ixs) ) = {0};
  fname_ixs = cellfun(@(ixs) ixs(end) + 1, fname_ixs);

  filenames = arrayfun( ...
    @(ix) filenames{ix}(fname_ixs(ix):end), ...
    1:length(fname_ixs), ...
    'UniformOutput', false);

  filenames = reshape(filenames, numel(filenames), 1);

  re_labels = re_patterns_table.Properties.VariableNames;
  n_patterns = numel(re_labels);

  function fname_matches = matches_in_filename(fname)
    fname_matches = cell(1, n_patterns);

    for px = 1:n_patterns
      p_label = re_labels{px};
      next_match = regexp(fname, char(re_patterns_table.(p_label) ), 'match');

      if isempty(next_match)
        fname_matches{px} = '';
      else
        fname_matches{px} = next_match{1};
      end
    end
  end

  fname_info = cellfun( ...
    @(fname) [matches_in_filename(fname), {fname}], ...
    filenames, ...
    'UniformOutput', false);

  fname_info = cell2table(cat(1, fname_info{:}) );
  fname_info.Properties.VariableNames = [re_labels, {'filename'}];
end

