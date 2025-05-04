function lfp_struct_mean = average_LFP(lfp_struct, chans_info)

  assert(isstruct(lfp_struct) );

  lfp_regions = fieldnames(lfp_struct);

  if ~exist('chans_info', 'var') || isempty(chans_info)
    chans_info.CHANNAMES = lfp_regions;
    chans_info.CHANACTIVE = true([numel(lfp_regions), 1]);

  else
    expected_chans_fields = [{'CHANNAMES'}, {'CHANACTIVE'}];

    % Make sure that the 'chans_info' struct has the correct fields.
    assert(isempty(setdiff(fieldnames(chans_info), expected_chans_fields) ));
    assert(isempty(setdiff(expected_chans_fields, fieldnames(chans_info) )));

    % Make sure the channels listed in the structure obtained from the
    % chans file match the channels and their order in the LFP structure.

    % The intersect() function finds the elements that occur in both lists,
    % and returns the indices needed to put both lists in order.
    [~, lfp_ixs, chans_ixs] = intersect( ...
      lfp_regions, ...
      chans_info.CHANNAMES);

    assert(numel(lfp_regions) == numel(chans_ixs) );

    % This finds the order of indices into `chans_info.CHANNAMES` such 
    % that the channels are ordered the same in `lfp_struct` but without
    % reordering the fields in `lfp_struct`.
    [~, sort_ixs] = sort(lfp_ixs);
    chans_ixs = chans_ixs(sort_ixs);

    chans_info.CHANNAMES = chans_info.CHANNAMES(chans_ixs);
    chans_info.CHANACTIVE = logical(chans_info.CHANACTIVE(chans_ixs) );

    assert(all(strcmp(lfp_regions, chans_info.CHANNAMES) ));
  end

  lfp_struct = double(cell2mat(struct2cell(lfp_struct) ));
  lfp_struct = lfp_struct(chans_info.CHANACTIVE,:);
  lfp_regions = lfp_regions(chans_info.CHANACTIVE);

  base_region_re_pattern = '(?i)^[_A-Z]*(?=((_[0-9]+$)|$))';
  base_regions = cellfun( ...
    @(chan) regexp(chan, base_region_re_pattern, 'match'), ...
    lfp_regions, ...
    'UniformOutput', true);

  [group_ixs, regions] = findgroups(base_regions);
  lfp_struct_mean = cell2struct( ...
    mat2cell( ...
      splitapply( ...
        @(group_mat) mean(group_mat, 1, 'omitnan'), ...
        lfp_struct, ...
        group_ixs), ...
      ones(1, numel(regions) )), ...
    regions);
end
