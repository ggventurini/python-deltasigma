function [snr,amp] = simulateQSNR(ntf,R,amp,f0,nlev,f,k)
%[snr,amp] = simulateQSNR(ntf,R,amp,f0=0,nlev=2,f=1/(4*R),k=13)
%Determine the SNR for a quadrature delta-sigma modulator using simulations.
%The modulator is described by a noise transfer function (ntf)
%and the number of quantizer levels (nlev).
%The ntf/stf may be given in ABCD form.
%The band of interest is defined by the oversampling ratio (R)
%and the center frequency (f0).
%The input signal is characterized by the amp vector and the f variable.
%amp defaults to [-120 -110...-20 -15 -10 -9 -8 ... 0]dB, where 0 dB means
%a full-scale (peak value = nlev-1) sine wave.
%f is the input frequency, normalized such that 1 -> fs;
%f is rounded to an FFT bin.
%
%Using sine waves located in FFT bins, the SNR is calculated as the ratio
%of the sine wave power to the power in all in-band bins other than those
%associated with the input tone. Due to spectral smearing, the input tone
%is not allowed to lie in bins 0 or 1. The length of the FFT is 2^k.
%
% Future versions may accommodate STFs.

% Handle the input arguments
if nargin<1
    error('Insufficient arguments');
end
parameters = {'ntf';'R';'amp';'f0';'nlev';'f';'k'};
defaults = [ NaN 64 NaN 0 2 NaN 13];
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults(i);'])
    end
end
if isnan(amp)
    amp = [-120:10:-20 -15 -10:0];
end
if isnan(f)
    f = f0 + 1/(4*R);
end

if abs(f-f0) > 1/(2*R)
    fprintf(1,'Warning: the input tone is out-of-band.\n');
end

N = 2^k;
if N < 8*R	% Require at least 8 bins to be "in-band"
    fprintf(1,'Warning: Increasing k to accommodate a large oversampling ratio.\n');
    k = ceil(log2(8*R))
    N = 2^k;
end
F = round(f*N);
if F<=1
    fprintf(1,'Warning: Increasing k to accommodate a low input frequency.\n');
    % Want f*N >= 1
    k = ceil(log2(1/f))
    N = 2^k;
    F = 2;
end

Ntransient = 100;
tone = (nlev-1) * exp(2i*pi*F/N*[-Ntransient:(N-1)]);
window = .5*(1 - cos(2*pi*(0:N-1)/N) );		%Hann window
f1 = max(round(N*(0.5+f0-0.5/R))+1,1);
inBandBins = [f1:round(N*(0.5+f0+0.5/R))+1];
F = F-f1+1+N/2;

snr = zeros(size(amp));
i=1;
for A = 10.^(amp/20);
    v = simulateQDSM(A*tone, ntf, nlev);
    hwfft = fftshift( fft(window.*v(1+Ntransient:N+Ntransient)) );
    snr(i) = calculateSNR(hwfft(inBandBins),F);
    i=i+1;
end
