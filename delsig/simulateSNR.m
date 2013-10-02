function [snr,amp] = simulateSNR(arg1,osr,amp,f0,nlev,f,k,quadrature)
%[snr,amp] = simulateSNR(ntf|ABCD|function,osr,amp,f0=0,nlev=2,f=1/(4*osr),k=13,quadrature=0)
%Determine the SNR for a delta-sigma modulator by using simulations.
%The modulator is described by a noise transfer function (ntf)
%and the number of quantizer levels (nlev).
%Alternatively, the first argument to simulateSNR may be an ABCD matrix or 
%the name of a function taking the input signal as its sole argument.
%The band of interest is defined by the oversampling ratio (osr)
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
% If ntf is complex, simulateQDSM (which is slow) is called.
% If ABCD is complex, simulateDSM is used with the real equivalent of ABCD
% in order to speed up simulations.

% Future versions may accommodate STFs.

% Handle the input arguments
if nargin<1
    error('Insufficient arguments');
end
parameters = {'arg1';'osr';'amp';'f0';'nlev';'f';'k';'quadrature'};
defaults = {NaN 64 NaN 0 2 NaN 13 0};
for arg_i=1:length(defaults)
    parameter = char(parameters(arg_i));
    if arg_i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
	 eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
	    eval([parameter '=defaults{arg_i};'])
    end
end
% Look at arg1 and decide if the system is quadrature
quadrature_ntf = 0;
if ischar(arg1)
    is_function = 1;
else
    is_function = 0;
    if isobject(arg1) % zpk object
        pz = [arg1.p{1} arg1.z{1}];
        for i=1:2
            if any( abs(imag(poly(pz(:,i)))) > 1e-4 )
                quadrature = 1;
        		quadrature_ntf = 1;
            end
        end
    else            % ABCD matrix
        if ~all(all(imag(arg1)==0))
            quadrature = 1;
        end
    end
end
if isnan(amp)
    amp = [-120:10:-20 -15 -10:0];
end
osr_mult = 2;
if f0~=0 & ~quadrature
    osr_mult = 2*osr_mult;
end
if isnan(f)
    f = f0 + 0.5/(osr*osr_mult);    % Halfway across the band
end
M = nlev-1;
if quadrature & ~quadrature_ntf	
    % Modify arg1 (ABCD) and nlev so that simulateDSM can be used
    nlev = [nlev; nlev];
    arg1 = mapQtoR(arg1);
end

if abs(f-f0) > 1/(osr*osr_mult)
    fprintf(1,'Warning: the input tone is out-of-band.\n');
end

N = 2^k;
if N < 8*2*osr	% Require at least 8 bins to be "in-band"
    fprintf(1,'Warning: Increasing k to accommodate a large oversampling ratio.\n');
    k = ceil(log2(8*2*osr))
    N = 2^k;
end
F = round(f*N);
if abs(F)<=1
    fprintf(1,'Warning: Increasing k to accommodate a low input frequency.\n');
    % Want f*N >= 1
    k = ceil(log2(1/f))
    N = 2^k;
    F = 2;
end

Ntransient = 100;
soft_start = 0.5*(1-cos(2*pi/Ntransient*[0:Ntransient/2-1])); 
if ~quadrature
    tone = M * sin(2*pi*F/N*[0:(N+Ntransient-1)]);
    tone(1:Ntransient/2) = tone(1:Ntransient/2) .* soft_start;
else
    tone = M * exp(2i*pi*F/N*[0:(N+Ntransient-1)]);
    tone(1:Ntransient/2) = tone(1:Ntransient/2) .* soft_start;
    if ~quadrature_ntf
    	tone = [real(tone); imag(tone)];
    end
end
window = .5*(1 - cos(2*pi*(0:N-1)/N) );		%Hann window
if f0==0
    % Exclude DC and its adjacent bin
    inBandBins = N/2+[3:round(N/(osr_mult*osr))];
    F = F-2;
else
    f1 = round(N*(f0-1/(osr_mult*osr)));
    inBandBins = N/2+[f1:round(N*(f0+1/(osr_mult*osr)))]; % Should exclude DC
    F = F-f1+1;
end

snr = zeros(size(amp));
i=1;
for A = 10.^(amp/20);
    if is_function
    	v = feval(arg1, A*tone);
    elseif quadrature_ntf
    	v = simulateQDSM(A*tone, arg1, nlev);
    else
    	v = simulateDSM(A*tone, arg1, nlev);
        if quadrature
    	    v = v(1,:) + 1i*v(2,:);
        end
    end
    hwfft = fftshift(fft(window.*v(1+Ntransient:N+Ntransient)));
    snr(i) = calculateSNR(hwfft(inBandBins),F);
    i=i+1;
end
