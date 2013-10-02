function y=dbm(v,R)
% dbm(v,R=50) = 10*log10(v^2/R*1000)  The equivalent in dBm of an rms voltage v
if isempty(v)
    return
end
if nargin<2
    R = 50;
end
y = -Inf*ones(size(v));
nonzero = v~=0;
y(nonzero) = 10*log10(abs(v(nonzero).^2)/R)+30;
