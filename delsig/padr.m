function y = padr(x, n, val)
%y = padr(x, n, val)
% Pad a matrix x on the right to length n with value val(0)
% The empty matrix is assumed to be have 1 empty row
if nargin < 3
    val = 0;
end

y = [ x val(ones(max(1,size(x,1)),n-size(x,2))) ];
