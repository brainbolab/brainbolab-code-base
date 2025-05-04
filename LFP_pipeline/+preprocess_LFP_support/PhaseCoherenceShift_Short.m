function [TheOffset] = PhaseCoherenceShift_Short( ...
  LFP_1, ...
  LFP_2, ...
  TheFreq, ...
  Fs)
% INPUT:
%   LFP_1: first LFP signal
%   LFP_2: second LFP signal
%   TheFreq: frequency of interest
%   Fs: sampling frequency
% OUTPUT:
%   TheOffset: optimized phase lag to maximize mean resultant length
% NOTE:
% If offset is positive, LFP2 leads LFP1
% if offset is negative, LFP1 leads LFP2

  % convert to double, if they weren't already
  LFP_1 = double(LFP_1);
  LFP_2 = double(LFP_2);

  % can change this if needed, but probably ok for resolution (ms)
  TimeStep = 2;
  % sampling stuff
  MaxOffset = Fs / 4;

  % create shift vector
  ShiftValues = (-MaxOffset):TimeStep:MaxOffset;

  % find out which signals to use
  AnalyticalSignal.LFP_1 = GetAnalyticalSignal(LFP_1);
  AnalyticalSignal.LFP_2 = GetAnalyticalSignal(LFP_2);

  LFP_1_Mask(1,1:length(LFP_1) ) = 0;
  LFP_2_Mask(1,1:length(LFP_2) ) = 0;

  LFP_1_Mask(MaxOffset+1:length(LFP_1_Mask)-MaxOffset) = ...
    AnalyticalSignal.LFP_1.Amplitude(MaxOffset+1:length(LFP_1_Mask)-MaxOffset) < 6000;

  LFP_2_Mask(MaxOffset+1:length(LFP_1_Mask)-MaxOffset) = ...
    AnalyticalSignal.LFP_2.Amplitude(MaxOffset+1:length(LFP_1_Mask)-MaxOffset) < 6000;

  % only use indices where both LFP_1 and LFP_2 are ok
  TheIdx = find(LFP_1_Mask & LFP_2_Mask);

  % build a butterworth filter with frequency band centered at the frequency
  [b, a] = butter(3, [(TheFreq-.5)/(Fs/2) (TheFreq+.5)/(Fs/2)]);

  % filter both signals
  LFP_1_Shift = filtfilt(b, a, LFP_1);
  LFP_2_Shift = filtfilt(b, a, LFP_2);

  % get the analytical signal for both signals
  AnalyticalSignal.LFP_1_Shift = GetAnalyticalSignal(LFP_1_Shift);
  AnalyticalSignal.LFP_2_Shift = GetAnalyticalSignal(LFP_2_Shift);

  % initialize variable for all of the shift magnitudes
  FreqShift = zeros(1, length(ShiftValues) );

  % loop through shifts
  for ShiftIdx = 1:length(ShiftValues)
	  ShiftVal = ShiftValues(ShiftIdx);
    % keep LFP 1 the same, shift LFP to the left
    
	  LFP_1_Seg = AnalyticalSignal.LFP_1_Shift.Phase(TheIdx);
	  LFP_2_Seg = AnalyticalSignal.LFP_2_Shift.Phase(TheIdx - ShiftVal); %501 to the end   gives amplitude and phase
    
    % take mrl
	  FreqShift(ShiftIdx) = circ_r(LFP_1_Seg - LFP_2_Seg);
  end %ShiftVal

  %% final result:
  % if shift is positive, LFP2 leads LFP1
  % if shift is negative, LFP1 leads LFP2
  ShiftValues = ShiftValues*1000 / Fs; 
  [FinalScore b] = max(FreqShift);
  TheOffset = ShiftValues(b);
end
