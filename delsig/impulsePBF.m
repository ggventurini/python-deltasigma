function [hc, t] = impulsePBF(C,x0,n,d)
%[hc, t] = impulsePBF(C,x0=0,n=20,d=0) Impulse response of a continuous-time polynomial-based filter (PBF)
%Input
% C     C(:,j) are the coefficients of segment j starting with the constant term
% x0    offset on x
% n     n points per segment
% d     derivative d
%
%Output
% hc    1x(n*size(C,2)) vector of the samples of the dth derivative of hc

% Argument checking and default-setting 
ArgumentsAndDefaults = {
  'C'   NaN 
  'x0'	0
  'n'	20
  'd'	0 
  };
for i = 1:size(ArgumentsAndDefaults,1)
    parameter = ArgumentsAndDefaults{i,1};
    if i>nargin || eval(['isempty(' parameter ') ']) || ...
      ( eval(['isnumeric(' parameter ') '])  &&  ...
        eval(['length(' parameter ') <= 1']) && ...
        eval(['isnan(' parameter ')']))
        if isnan(ArgumentsAndDefaults{i,2})
            error('%s: Argument %d (%s) is required.',mfilename, i, parameter )
        else
            eval([parameter '= ArgumentsAndDefaults{i,2};'])
        end
    end
end

N = size(C,1);
if nargout > 1
    t = linspace(0,N-1/n,n*N);
end

M = size(C,2)-1;
if M<d
    hc = zeros(size(t));
else
    for i = 1:d
        C = C(:,2:end).*repmat(1:M,N,1);
        M = M-1;
    end
    x = linspace(0,1-1/n,n) + x0;
    X = repmat(x,M+1,1).^repmat((0:M)',1,n);
    hc = (C*X)';
    hc = hc(:).';
end
