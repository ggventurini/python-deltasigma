function y = circ_smooth(x,n)
nx = length(x);
w = ds_hann(n)/(n/2);
xw = conv(x,w);
y = circshift([ xw(n:nx) xw(1:n-1)+xw(nx+1:end)],[0 n/2-1]);
