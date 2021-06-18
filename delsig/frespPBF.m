function Hc = frespPBF(f,C,x0)
%Hc = frespPBF(f,C,x0=-0.5)	
%Compute the frequency response of a polynomial-based filter.
% Input
%  f    A vector of frequencies.
%  C    Nx(M+1) matrix containing the coefficients of the polynomial pieces.
%       p_i(x) = C(i,1) + C(i,2)*x + C(i,3)*x^2 ... + C(i,M+1)*x^M
%  x0	Offset on the polynomial argument, i.e. x = mu + x0 where mu is in [0,1].
%
% Output
%  Hc   Hc(f)

% Handle the input arguments
% Argument checking and default-setting
ArgumentsAndDefaults = {
    'f'     NaN
    'C'     NaN
    'x0'    -0.5
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
f = f(:)'; % Make f a row vector

[N Mp1] = size(C);
M = Mp1 - 1;
nf = length(f);
nz = f~=0;
w = 2*pi*f;
% F(n+1,j) is the Fourier transform of one x^m segment evaluated at w_j.
F = zeros(M+1,nf);
% Handle f=0
mp1 = 1:Mp1;
if any(~nz)
    F(:,~nz) = ((x0+1).^mp1 - x0.^mp1) ./ mp1;
end
% Handle nonzero frequencies
wnz = w(nz);
F(1,nz) = sinc(f(nz)).*exp(-1i*(x0+0.5)*wnz);
exp1 = exp(-1i*(x0+1)*wnz);
exp0 = exp(-1i*(x0)*wnz);
for m=1:M
    F(m+1,nz) = 1i./wnz.*( (x0+1)^m*exp1 - x0^m*exp0 - m*F(m,nz) );
end
p = exp(-1i*w);
ph = ones(size(w));
Hc = 0;
for n = 1:N
    Hc = Hc + C(n,:)*F.*ph;
    ph = ph .* p;
end

