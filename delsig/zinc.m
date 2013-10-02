function mag = zinc(f,m,n)
% mag = zinc(f,m=64,n=1)	Calculate the magnitude response
% of a cascade of n mth-order comb filters at frequencies f.

if nargin<3
    n = 1;
    if nargin<2
	m = 64;
    end
end

mag = abs( sinc(m*f) ./ sinc(f) ).^n;
