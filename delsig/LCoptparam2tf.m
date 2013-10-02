function [H,L0,ABCD,param] = LCoptparam2tf(x,param)
% LCoptparam2tf(...)	Convert optimization parameters to transfer functions.

n = param.n;
param.gw = ones(1,n-1);
switch param.form
  case 'FB'
    param.gu = [1 zeros(1,n)];
    param.gv = x;
    param.gx = [zeros(1,n-1) 1];
    param.rx = zeros(1,n);
  case 'FF'
    param.gu = [1 zeros(1,n-1) 1];
    param.gv = [1 zeros(1,n)];
    param.gx = x;
    param.rx = zeros(1,n);
  case 'FFR'
    param.gu = [1 zeros(1,n-1) 1];
    param.gv = [x(1) zeros(1,n)];
    param.gx = [zeros(1,n-1) 1];
    param.rx = [0 x(2:n)];
  case 'GEN'
    param.gv = [x(1:n)];
    param.rx = [0 x(n+1:2*n-1)];
  otherwise
    fprintf(2,'%s: Error. form=%s is not supported.\n',mfilename,form);
    return
end

[H L0 ABCD] = LCparam2tf(param);
