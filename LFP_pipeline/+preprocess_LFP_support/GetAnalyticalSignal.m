function AnalyticalSignal = GetAnalyticalSignal(TheLFP)
  % transform the LFP using Hilbert
  HilbertLFP = hilbert(TheLFP);
  % get the phase array from the tangent of the imaginary and real portions
  HilbertPhaseArray = atan2(imag(HilbertLFP), real(HilbertLFP) );
  % get amplitude using pythagorean theorem
  HilbertAmplitudeArray = sqrt((real(HilbertLFP).^2 + imag(HilbertLFP).^2)');

  % save information and return
  AnalyticalSignal.Phase = HilbertPhaseArray;
  AnalyticalSignal.Amplitude = HilbertAmplitudeArray;
end