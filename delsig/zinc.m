function mag = zinc(f,n,m)
% mag = zinc(f,n=64,m=1)
% Calculate the magnitude response of a cascade of m sinc_n filters at frequencies f.
if nargin<3
    m = 1;
    if nargin<2
	n = 64;
    end
end

mag = abs( sinc(n*f) ./ sinc(f) ).^m;
