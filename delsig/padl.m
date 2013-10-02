function y = padl(x, n, val)
%y = padl(x, n, val)
% Pad a matrix x on the left to length n with value val(0)
% The empty matrix is assumed to be have 1 empty row
if nargin < 3
    val = 0;
end

y = [ val(ones(max(1,size(x,1)),n-size(x,2))) x ];
