function [peak_snr,peak_amp] = peakSNR(snr,amp) 
% [peak_snr,peak_amp] = peakSNR(snr,amp)	Find the snr peak by fitting
% a smooth curve to the top end of the SNR vs input amplitude data.
% Both amp and snr are expressed in dB.

% Delete garbage data
i = Inf == abs(snr); snr(i) = []; amp(i) = [];
i = isnan(snr); snr(i) = []; amp(i) = [];

max_snr = max(snr);
i = find(snr>max_snr-10);
min_i = min(i); max_i = max(i);
j = find( snr(min_i:max_i) < max_snr-15 );
if j
	max_i = min_i + min(j) -2;
	i = min_i:max_i;
end
snr = 10.^(snr(i)/20);
amp = 10.^(amp(i)/20);
n = length(i);

c = max(amp)*1.05;
% fit y = ax + b/(x-c) to the data
A = [amp; 1./(amp-c)];
ab = snr/A; 

peak_amp = c- sqrt(ab(2)/ab(1));
peak_snr = ab*[peak_amp; 1/(peak_amp-c)];
peak_snr = dbv(peak_snr);
peak_amp = dbv(peak_amp);

return
% For debugging
amp = linspace(min(amp),max(amp));
pred = ab*[amp; 1./(amp-1)];
plot(amp,pred,'g')

