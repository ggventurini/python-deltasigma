function y = delay(x,n)
% y = delay(x,n=1) Delay signal x by n samples
if nargin < 2
    n = 1;
end
nx = length(x);
if nx <= n
    y = zeros(size(x));
else
    if size(x,1) > size(x,2)
	y = [zeros(n,1); x(1:nx-n)];
    else
	y = [zeros(1,n) x(1:nx-n)];
    end
end
