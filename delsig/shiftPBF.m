function Cp = shiftPBF(C,x0)
% Cp = shiftPBF(C,x0) 
% Shift the polynomial argument of a polynomial-based filter by x0,
% i.e. x' = x + x0

[N, Mp1] = size(C);
Cp = zeros(N,Mp1);
M = Mp1 - 1;
m = 0:M;
x = (x0.^m)';
for i=1:size(C,2)
    Cp(:,i) = C*x;
    C = C .* repmat(m,N,1);
    C = C(:,2:end);
    m(end) = [];
    x(end) = [];
end
