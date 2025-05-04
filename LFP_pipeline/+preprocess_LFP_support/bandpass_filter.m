function results = bandpass_filter(signal, fs, high_freq, low_freq)
% bandpass filter
  high_pass = 2 * (high_freq / fs);
  low_pass  = 2 * (low_freq / fs);
  [butter_b, butter_a] = butter(2, [low_pass, high_pass]);
    
  results = filtfilt(butter_b, butter_a, signal);
end