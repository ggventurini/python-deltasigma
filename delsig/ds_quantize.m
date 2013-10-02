function v = ds_quantize(y,n)
%v = ds_quantize(y,n=2)
%Quantize y to 
% an odd integer in [-n+1, n-1], if n is even, or
% an even integer in [-n+1, n-1], if n is odd.
%This definition gives the same step height for both mid-rise
%and mid-tread quantizers.
%n can be a column vector which specifies how to quantize the rows of y
if nargin<2
    n=2;
end
if length(n)==1
    n = n*ones(size(y));
elseif size(n,1)==size(y,1)
    n = n*ones(1,size(y,2));
end

i = rem(n,2)==0;	
v = zeros(size(y));
v(i) = 2*floor(0.5*y(i))+1;	% mid-rise quantizer
v(~i) = 2*floor(0.5*(y(~i)+1));	% mid-tread quantizer

% Limit the output
L = n-1;
for m=[-1 1]
    i = m*v>L;
    if any(i(:))
	v(i) = m*L(i);
    end
end
