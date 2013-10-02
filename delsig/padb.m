function y = padb(x, n, val)
%y = padb(x, n, val)
% Pad a matrix x on the bottom to length n with value val(0)
% The empty matrix is assumed to be have 1 empty column
if nargin < 3
    val = 0;
end

y = [ x ; val( ones( n-size(x,1), max(1,size(x,2)) ) ) ];
