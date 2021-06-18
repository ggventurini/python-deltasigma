function plotSpectrum(X,fin,fmt,nbin,n)
% plotSpectrum(X,fin,fmt='-',nbin=8,n=3) Plot a smoothed spectrum
if nargin<5
    n = 3;
    if nargin<4
        nbin = 8;
        if nargin<3
            fmt = '-';
        end
    end
end
i = fmt=='.' | ( fmt>='0' & fmt <='9' );
if any(i)
    i = find(~i); i = i(end); % i = index of last character before the numeric suffix
    lw = str2num(fmt(i+1:end));
    fmt = fmt(1:i);
else
    lw = 0.5;
end
[f p] = logsmooth(X,fin,nbin,n);
semilogx(f,p,fmt,'Linewidth',lw);
