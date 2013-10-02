function y = rms(x,no_dc)
% y = rms(x,no_dc(0))
if nargin<2
    no_dc=0;
end
if no_dc
    x = x - mean(x);
end
y = norm(x)/sqrt(length(x));

