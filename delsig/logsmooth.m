function [f,p] = logsmooth(X,inBin,nbin,n)
% function [f,p] = logsmooth(X,inBin,nbin=8,n=3)
% Smooth the fft, X, and convert it to dB.
% Use nbin bins from 0 to 3*inBin, 
% thereafter increase bin sizes by a factor of 1.1, staying less than 2^10.
% For the n sets of bins inBin+[0:2], 2*inBin+[0:2], ... 
% n*inBin+[0:2], don't do averaging. This way, the noise BW
% and the scaling of the tone and its harmonics are unchanged.
% Unfortunately, harmonics above the nth appear smaller than they 
% really are because their energy is averaged over many bins.
parameters = {'X','inBin','nbin','n'};
defaults = { [] [] 8 3 };
for arg_i=1:length(defaults)
    parameter = char(parameters(arg_i));
    if arg_i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{arg_i};'])
    end
end
N=length(X); N2=floor(N/2);
f1 = rem(inBin-1,nbin)+1;
startbin = [f1:nbin:inBin-1 inBin:inBin+2];
for i=1:n
	startbin = [startbin startbin(end)+1:nbin:i*inBin-1 i*inBin+[0:2] ];
end
m = startbin(length(startbin))+nbin;
while m < N2
   startbin = [startbin m];
   nbin = min(nbin*1.1,2^10);
   m = round(m+nbin);
end
stopbin = [startbin(2:length(startbin))-1 N2];
f = ((startbin+stopbin)/2 -1)/N;
p = zeros(size(f));
for i=1:length(f)
    p(i) = dbp(norm(X(startbin(i):stopbin(i)))^2 / ...
              (stopbin(i)-startbin(i)+1));
end
