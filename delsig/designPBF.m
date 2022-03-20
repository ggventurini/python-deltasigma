function [C, e, x0] = designPBF(N,M,pb,pbr,sbr,ncd,np,ns,fmax)
% [C, e, x0] = designPBF(N=10,M=5,pb=0.25,pbr=0.1,sbr=-100,ncd=0,np=100,ns=1000,fmax=5)
% Design a symmetric polynomial-based filter (PBF) according to Hunter's method [1].
% Inputs
%  N    Number of polynomial pieces.
%  M	Order of the polynomial pieces.
%  pb   Passband width. Relative the input sample rate,
%       the passband is [0,pb] and the stopband is [1-pb,Inf].
%  pbr  Passband ripple in dB.
%  sbr  Stobpand ripple in dB.
%  ncd	Number of continuous derivatives. -1 if the impulse response itself is not continuous
%  np   Number of points in the passband.
%  ns   Number of points in the stopband.
%  fmax Maximum frequency checked in the stopband.
%
% Output
%  C    Nx(M+1) matrix containing the coefficients of the polynomial pieces.
%       p_i(x) = C(i,1) + C(i,2)*x + C(i,3)*x^2 ... + C(i,M+1)*x^M
%  e    The maximum weighted error. e<=1 implies the specs were met.
%  x0	Offset on the polynomial argument, i.e. x = mu + x0 where mu is in [0,1].
%       In this implementation, x0 is always -0.5.
%
% [1] M. T. Hunter, "Design of polynomial-based filters for continuously
% variable sample rate conversion with applications in synthetic
% instrumentation and software defined radio," Ph.D. thesis,
% University of Florida, 2008.
%
%Future Enhancements
% Allow N to be odd.
% Support arbitrary passband functions and weights

% Argument checking and default-setting
ArgumentsAndDefaults = {
    'N'     10
    'M'     5
    'pb'    0.25
    'pbr'   0.1     % dB
    'sbr'   -100    % dB
    'ncd'   0
    'np'    100
    'ns'    1000
    'fmax'  5
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
odd = mod(N,2);
N2 = (N-odd)/2;
x0 = -0.5;
if pb >= 0.5
    error('pb, the passband width, must be less than 0.5.');
end

% Frequencies at which Hc(f) is calculated.
w = 2*pi*[linspace(0,pb,np) linspace(1-pb,fmax,ns)];

pbwt = 1/(undbv(pbr)-1);
sbwt = 1/undbv(sbr);
invWt = [repmat(1/pbwt,np,1); repmat(1/sbwt,ns,1)];
% Construct matrices for linprog
f = zeros( (N2+odd)*(M+1)+1, 1 ); f(end) = 1;
FT = computeFT(w,N,M);
A = [ FT -invWt; -FT -invWt ];
d = [ones(np,1); zeros(ns,1)];
b = [d; -d];
Aeq = equality_constraints(N,M,ncd);
beq = zeros(size(Aeq,1),1);
x = initialGuess(N,M);
% Solve linear program
options.Display = 'off';
[x, e, exitflag, output] = linprog(f,A,b,Aeq,beq,[],[],x,options); 
if exitflag <= 0
    fprintf(1,'Sorry. linprog returned an error:\n %s\n',output.message);
end
% Construct the C matrix
if odd % Inconsistent! My odd and even formulations populate Cb differently
    Cb = reshape(x(1:end-1),M+1,N2+odd)';
else
    Cb = reshape(x(1:end-1),N2,M+1);
end
Ct = flipud(Cb((1+odd):end,:)).*repmat((-1).^(0:M),N2,1);
C = [Ct; Cb];
end

% Compute a matrix which gives the Fourier Transform at the specified
% frequencies from the polynomial coefficient vector
function FT = computeFT(w,N,M)
odd = mod(N,2);
N2 = (N-odd)/2;
nf = length(w);
nz = w~=0;
% F(n+1,j) is the Fourier transform of one x^m segment evaluated at w_j.
F = zeros(M+1,nf);
% Handle zero frequencies separately
m = 1:M+1;
if any(~nz)
    F(:,~nz) = (0.5.^m - (-0.5).^m) ./ m;
end
% Handle nonzero frequencies
wnz = w(nz);
F(1,nz) = sinc(wnz/(2*pi));
for m=1:M
    if rem(m,2)==1
        F(m+1,nz) = 1i./wnz.*( 0.5^m*2*cos(wnz/2) - m*F(m,nz));
    else
        F(m+1,nz) = 0.5^m*F(1,nz) + m./(1i*wnz).*F(m,nz);
    end
end
if odd % Inconsistent! My odd and even formulations populate x differently
    % p contains the phase factors for each frequency.
    % p = exp((-i*w_j)*n), n=0:N2
    F = F.';
    FT = zeros(nf,(N2+odd)*(M+1));
    m = 1:M+1;
    FT(:,m) = real(F);
    p1 = exp(-1i*w');
    p = p1;
    for i = 1:N2
        m = m + M+1;
        FT(:,m) = 2*real(repmat(p,1,M+1).*F);
        if i<N2
            p = p.*p1;
        end
    end
else
    % P contains the phase factors for each (segment,frequency).
    % P(n,j) = exp((-i*w_j)*(n-0.5))
    % This is more compact, but slower:
    % P = exp(repmat((0.5:(N-1)/2)',1,nf) .* repmat(-1i*w,N/2,1));
    P = zeros(N2,nf);
    P(1,:) = exp(-0.5i*w);
    p = P(1,:).^2;
    for n = 2:N2
        P(n,:) = P(n-1,:) .* p;
    end
    % See DSToolbox notebook page xxx
    FT = 2 * real( repmat(P,M+1,1) .* reshape( repmat(F(:).',N2,1), N2*(M+1), nf ) )';
end
end

function Aeq = equality_constraints(N,M,ncd)
odd = mod(N,2);
N2 = (N-odd)/2;
Mp1 = M+1;
if odd
    Aeq = zeros(ceil(M/2)+(N2+1)*(ncd+1), (N2+odd)*Mp1+1);
    % Odd powers of x in C0 must be zero
    for i = 1:ceil(M/2)
        Aeq(i,2*i) = 1;
    end
    % Continuity constraints
    m = 0:M;
    a1 = 0.5.^m;
    a2 = -(-0.5).^m;
    for d = 0:ncd
        for n = 0:N2
            i = i+1;
            Aeq(i,n*Mp1+(1:Mp1)) = a1;
            if n<N2
                Aeq(i,(n+1)*Mp1+(1:Mp1)) = a2;
            end
        end
        if d < ncd
            a1 = a1.*m*2;
            a2 = a2.*m*(-2);
            m = m - 1;
        end

    end
else
    Aeq = zeros(N2*(ncd+1)+floor(ncd/2), N2*(M+1)+1);
    if ncd >= 0
        a1 = 0.5.^(0:M);
        a2 = -(-0.5).^(0:M);
        i = 1:N2:M*N2+1;
        for n = 1:N2
            Aeq(n,i) = a1;
            i = i + 1;
            if n<N2
                Aeq(n,i) = a2;
            end
        end
        nn = N2+1;
        m = 0:M;
        for d = 1:ncd
            a1 = a1.*m*2;
            a2 = a2.*m*(-2);
            i = 1:N2:M*N2+1;
            if rem(d,2)==1 %odd-order derivatives have extra constraints at mu=0
                Aeq(nn,i) = a2;
                nn = nn + 1;
            end
            for n = 1:N2
                Aeq(nn,i) = a1;
                i = i + 1;
                if n<N2
                    Aeq(nn,i) = a2;
                end
                nn = nn + 1;
            end
            m = m - 1;
        end
    end
end
end

% Approximate sinc(n) with N polynomial pieces
function x = initialGuess(N,M)
odd = mod(N,2);
N2 = (N-odd)/2;
x = zeros((N2+odd)*(M+1)+1,1);
x(end) = 1;
Cb = zeros(N2+odd,M+1);
i = 1:N2+odd;
xi = i-(1+odd)/2;
Cb(i,1) = sinc(xi);
%Use finite differences because I am too lazy to code the derivatives
if M>0
    di = 0.1;
    Cb(i,2) = (sinc(xi+di/2) - sinc(xi-di/2)) / di;
    if M>1
        Cb(i,3) = (sinc(xi+di) - 2*sinc(xi) + sinc(xi-di)) / di^2/2;
    end
end
if odd
    x(1:end-1) = reshape(Cb',(N2+odd)*(M+1),1);
else
    x(1:end-1) = Cb(:);
end
end

