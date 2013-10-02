function snr = calculateSNR(hwfft,f,nsig)
% snr = calculateSNR(hwfft,f,nsig=1) Estimate the signal-to-noise ratio, 
% given the in-band bins of a (Hann-windowed) fft and 
% the location of the input signal (f>0).
% For nsig=1, the input tone is contained in hwfft(f:f+2); 
%  this range is appropriate for a Hann-windowed fft.
% Each increment in nsig adds a bin to either side.
% The SNR is expressed in dB.
if nargin<3
    nsig = 1;
end
signalBins = [f-nsig+1:f+nsig+1];
signalBins = signalBins(signalBins>0);
signalBins = signalBins(signalBins<=length(hwfft));
s = norm(hwfft(signalBins));		% *4/(N*sqrt(3)) for true rms value;
noiseBins = 1:length(hwfft);
noiseBins(signalBins) = [];
n = norm(hwfft(noiseBins));
if n==0
    snr = Inf;
else
    snr = dbv(s/n);
end
