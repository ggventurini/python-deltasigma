function [ABCDq, ABCDp] = mapR2Q( ABCDr )
% [ABCDq ABCDp] = mapR2Q( ABCDr )	Map a real ABCD to a quadrature ABCD
% ABCDr has its states paired (real,imaginary)
%
% ABCDq is the quadrature (complex) version of ABCDr
% ABCDp is the mirror-image system matrix 
% (ABCDp is zero if ABCDr has no quadrature errors)

ABCD11 = ABCDr(1:2:end,1:2:end);
ABCD12 = ABCDr(1:2:end,2:2:end);
ABCD21 = ABCDr(2:2:end,1:2:end);
ABCD22 = ABCDr(2:2:end,2:2:end);

ABCDq = 0.5*(ABCD11+ABCD22) + 0.5i*(ABCD21-ABCD12);
ABCDp = 0.5*(ABCD11-ABCD22) + 0.5i*(ABCD21+ABCD12);

