function plotSpectrum(X,fin,fmt)
% plotSpectrum(X,fin,fmt) Plot a smoothed spectrum
if nargin<3
    fmt='-';
end
[f p] = logsmooth(X,fin);
semilogx(f,p,fmt);
