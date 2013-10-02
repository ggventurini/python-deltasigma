function y = sinc_decimate(x,m,r)
% y = sinc_decimate(x,m,r)	Decimate x by m-th order sinc of length r.
x = x(:)';
for i=1:m
	x = cumsum(x);
	x = [x(1:r) x(r+1:end)-x(1:end-r)]/r;
end
y = x(r:r:end);
