function bilogplot(V,f0,fbin,x,y,fmt)
% bilogplot(V,f0,fbin,x,y)	Plot two side-by-side spectra
% V    Hann-windowed FFT
% f0   bin number of center frequency
% fbin bin number of test tone
% x    xmin, xmax_left, xmax_right
% y    ymin, ymax, dy, y_skip

Xl = V(f0:-1:1);
Xr = V(f0:end);
N = length(V)-1;
fbin = abs(fbin-f0);
[f,p] = logsmooth2(Xl,fbin);
f = (f0-1)/N*f;
subplot('position',[0.05 .1 .45 .8])
semilogx(f,p,fmt);
hold on;
set(gca,'Xdir','reverse')
axis([x(1) x(2) y(1) y(2)]);
grid on;
ytix = y(1):y(3):y(2);
set(gca,'YTick', ytix);
set(gca,'YTickLabel', axisLabels(ytix,y(4)));

[f,p] = logsmooth2(Xr,fbin);
f = (N-f0+1)/N*f;
subplot('position',[0.5 .1 .45 .8])
semilogx(f,p,fmt);
hold on;
axis([x(1) x(3) y(1) y(2)]);
grid on;
ytix = y(1):y(3):y(2);
set(gca,'YTick', ytix);
set(gca,'YAxisLocation', 'right');
set(gca,'YTickLabel', axisLabels(ytix,y(4)));

return


function [f,p] = logsmooth2(X,inBin,nbin)
% function [f,p] = logsmooth2(X,inBin,nbin)
% Smooth the fft, X, and convert it to dB.
% Use nbin(8) bins from 0 to 3*inBin, 
% thereafter increase bin sizes by a factor of 1.1, staying less than 2^10.
% For the 3 sets of bins inBin+[0:2], 2*inBin+[0:2], and 
% 3*inBin+[0:2], don't do averaging. This way, the noise BW
% and the scaling of the tone and its harmonics are unchanged.
% Unfortunately, harmonics above the third appear smaller than they 
% really are because their energy is averaged over several bins.
if nargin<3
    nbin = 8;
end
N=length(X);
n = nbin;
f1 = rem(inBin-1,n)+1;
startbin = [f1:n:inBin-1 inBin:inBin+2 inBin+3:n:2*inBin-1 ...
	 2*inBin+[0:2] 2*inBin+3:n:3*inBin-1 3*inBin+[0:2] ];
m = startbin(length(startbin))+n;
while m < N
   startbin = [startbin m];
   n = min(n*1.1,2^10);
   m = round(m+n);
end
stopbin = [startbin(2:length(startbin))-1 N];
f = ((startbin+stopbin)/2 -1)/N;
p = zeros(size(f));
for i=1:length(f)
    p(i) = dbp(norm(X(startbin(i):stopbin(i)))^2 / ...
              (stopbin(i)-startbin(i)+1));
end
