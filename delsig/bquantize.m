function y = bquantize(x,nsd,abstol,reltol)
%y = bquantize(x,nsd=3,abstol=eps,reltol=10*eps)
% Bidirectionally quantize a n by 1 vector x to nsd signed digits, 
% Terminate early if the error is less than the specified tolerances.
% y is a structure array with the following fields:
% y(i).val is the quantized value in floating-point form
% y(i).csd is a 2-by-nsd (or less) matrix containing
% the powers of two (first row) and their signs (second row).
% See also bunquantize.m.

% Handle the input arguments
parameters = {'x' 'nsd' 'abstol' 'reltol'};
defaults = { NaN 3 eps 10*eps };
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end

n = length(x);
q = zeros(2*n,nsd);
y = struct('val',cell(1,n),'csd',cell(1,n));
offset = -log2(0.75);

for i = 1:n
    xp = x(i);
    y(i).val = 0;
    for j = 1:nsd
	error = abs(y(i).val-x(i));
	if error <= abstol | error <= abs(x(i))*reltol
	    break;
	end
	p = floor( log2(abs(xp)) + offset );
	p2 = pow2( p );
	sx = sign(xp);
	xp = xp- sx*p2;
	y(i).val = y(i).val + sx*p2;
	y(i).csd(1:2,j) = [p;sx]; 
    end
end
