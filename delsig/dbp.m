function y=dbp(x)
% dbp(x) = 10*log10(x): the dB equivalent of the power x
y = -Inf*ones(size(x));
if isempty(x)
    return
end
nonzero = x~=0;
y(nonzero) = 10*log10(abs(x(nonzero)));
