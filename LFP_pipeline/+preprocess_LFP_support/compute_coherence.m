function [CXY, freq, region_labels, time] = compute_coherence( ...
  lfp_struct, ...
  in_region, ...
  out_regions, ...
  win_len, ...
  win_step, ...
  fs)
% Function call:
%   compute_coherence(lfp_struct, in_region, out_regions, win_len, win_step, fs)
%
% Description:
%   Compute the magnitude-squared coherence between the region indicated by 'in_region' and each region listed
%   in the 'out_regions' array whose data is stored in 'lfp_struct'. The 'win_len' and 'win_step' variables
%   are used to indicate window length (in seconds) and the amount of overlap (in seconds) between subsequent
%   time windows. If you do not want windows to overlap, set the 'win_len' and 'win_step' variables to be
%   the same.
%
% Arguments:
%    lfp_struct: MATLAB struct whose fields contain the LFP recordings for different brain regions.
%
%     in_region: Label of one of the regions corresponding to a field of 'lfp_struct' which will be paired
%                with all of the regions in 'out_regions' in order to compute coherence.
%
%   out_regions: Cell array containing the labels of the regions corresponding to fields of 'lfp_struct'
%                whose LFP recordings will be paired with the region indicated by 'in_region' in order to
%                compute coherence.
%
%       win_len: Integer indicating how many seconds of data to use in each window. This should be at least
%                as large as 'win_step'.
%
%      win_step: Integer indicating the step size (in seconds) between starting times of subsequent windows.
%                This should be less than or equal to 'win_len'.
%
%            fs: The sampling frequency of the LFP data.
%
% Example:
%   >>> lfp_struct = load("./LFP_data/Mouse1_250101_LFP.mat");
%   >>> in_region = 'Acc';
%   >>> out_regions = [{'LPBN'}, {'MD_thal'}, {'Po'}, {'RPBN'}, {'VPM'}];
%   >>> win_len = 1;
%   >>> win_step = 1;
%   >>> fs = 1000;
%   >>> [CXY, freq, region_labels, time] = compute_coherence(lfp_struct, in_region, out_regions, win_len, win_step, fs);

  % Check to make sure 'win_len', 'win_step' and 'fs' are greater than zero.
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

  if ischar(out_regions)
    out_regions = {out_regions};
  end

  out_regions = reshape(out_regions, 1, []);
  regions = reshape(fieldnames(lfp_struct), 1,[]);

  assert(ischar(in_region) );
  assert(all(ismember([in_region, out_regions], regions) ));

  ix_in = find(strcmp(regions, in_region) );
  ixs_out = find(ismember(regions, out_regions) );

  out_regions = regions(ixs_out);
  region_labels = cellfun( ...
    @(out_region) join([{in_region}, {out_region}], '_x_'), ...
    out_regions)';

  lfp_mat = double(cell2mat(struct2cell(lfp_struct) ))';
  lfp_mat = lfp_mat(:,[ix_in, ixs_out]);
  n_data_pts = size(lfp_mat, 1);
  n_blocks = floor(n_data_pts / fs);
  n_windows = floor((n_blocks - win_len) / win_step) + 1;
  nfft = win_len * fs;

  win_ixs_table = nan(2, n_windows);
  win_ixs_table(1,:) = (0:win_step:(n_blocks - win_len) ) * fs + 1;
  win_ixs_table(2,:) = (win_len:win_step:n_blocks) * fs;
  time = reshape(floor(win_ixs_table(1,:) / fs) + 1, [], 1);

  assert(all((win_ixs_table(2,:) - win_ixs_table(1,:) + 1) == nfft) );

  CXY = nan(floor(nfft / 2) + 1, length(out_regions), n_windows);
  freq = [];

  for wx = 1:n_windows
    ax = win_ixs_table(1,wx);
    bx = win_ixs_table(2,wx);
    [CXY_ix, freq_ix] = mscohere(lfp_mat(ax:bx,1), lfp_mat(ax:bx,2:end), [], [], nfft, fs, 'mimo');
    CXY(:,:,wx) = CXY_ix;

    if isempty(freq)
      freq = reshape(freq_ix, [], 1);
    end
  end
end
