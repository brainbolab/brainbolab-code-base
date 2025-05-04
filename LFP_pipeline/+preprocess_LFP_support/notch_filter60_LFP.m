function lfp_struct = notch_filter60_LFP(lfp_struct, fs)

  lfp_regions = fieldnames(lfp_struct);
  lfp_struct = double(cell2mat(struct2cell(lfp_struct) ))';

  % Notch filter 60Hz signal
  filt_freq = 2 * (60 / fs); % Filter out 60Hz Signal
  filt_bw = filt_freq / 60; % 1Hz Bandwidth
  [filt_b, filt_a] = iirnotch(filt_freq, filt_bw);
  lfp_struct = filtfilt(filt_b, filt_a, lfp_struct)';
  lfp_struct = cell2struct( ...
    mat2cell(lfp_struct, repelem(1, length(lfp_regions) )), ...
    lfp_regions);

end
