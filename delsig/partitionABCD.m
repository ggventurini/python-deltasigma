function [A, B, C, D] = partitionABCD(ABCD, m);
% function [A B C D] = partitionABCD(ABCD, m); Partition ABCD into
% A, B, C, D for an m-input state-space system.
if nargin<2
    n = min(size(ABCD))-1;
    m = size(ABCD,2)-n;
else
    n = size(ABCD,2)-m;
end
r = size(ABCD,1)-n;

A = ABCD(1:n, 1:n);
B = ABCD(1:n, n+1:n+m);
C = ABCD(n+1:n+r, 1:n);
D = ABCD(n+1:n+r, n+1:n+m);
