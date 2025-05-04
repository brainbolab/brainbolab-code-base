function filenames = ls_filenames_w_pattern(dir_path, re_pattern)
  dir_path = char(dir_path);
  re_pattern = char(re_pattern);

  file_sep = '/';
  if ispc
    file_sep = '\';
  end

  if (dir_path(end) ~= file_sep)
    dir_path = [dir_path, file_sep];
  end

  filenames_all = split(ls(dir_path) );
  fname_ixs = cellfun(@(fname) ~isempty(fname), filenames_all);
  filenames_all = filenames_all(fname_ixs);

  match_ixs = true([numel(filenames_all), 1]);

  if ~isempty(re_pattern)
    match_ixs = cellfun( ...
      @(mtch) ~isempty(mtch), ...
      regexp(filenames_all, re_pattern) );
  end

  filenames = filenames_all(match_ixs);
end