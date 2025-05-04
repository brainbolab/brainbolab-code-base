function phase_offsets = compute_directionality( ...
  lfp_struct, ...
  region_1, ...
  region_2, ...
  freq_range, ...
  fs)
% Function call:
%   compute_directionality(lfp_struct, region_1, region_2, freq_range, fs)
%
% Description:
%   X
%
% Arguments:
%   lfp_struct: MATLAB struct whose fields contain the LFP recordings for different brain regions.
%
%     region_1: X
%
%     region_2: X
%
%   freq_range: X
%
%            fs: The sampling frequency of the LFP data.
%
% Example:
%   >>> lfp_struct = load("./LFP_data/Mouse1_250101_LFP.mat");
%   >>> region_1 = 'Acc';
%   >>> region_2 = 'LPBN';
%   >>> freq_range = [1, 50];
%   >>> fs = 1000;
%   >>> phase_offsets = compute_directionality(lfp_struct, region_1, region_2, freq_range, fs);

  assert(numel(freq_range) == 2);
  assert(mod(freq_range(1), 1) == 0);
  assert(mod(freq_range(2), 1) == 0);

  region_1 = char(region_1);
  region_2 = char(region_2);

  freq_range = freq_range(1):freq_range(2);
  n_freq = length(freq_range);

  phase_offsets = NaN(n_freq, [], 1);

  parfor(fx = 1:n_freq, 40)
    freq_fx = freq_range(fx);
    phase_offsets(fx) = PhaseCoherenceShift_Short( ...
      reshape(lfp_struct.(region_1), [], 1), ...
      reshape(lfp_struct.(region_2), [], 1), ...
      double(freq_fx), ...
      fs);
  end
end
