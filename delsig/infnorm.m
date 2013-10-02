function [Hinf,wmax] = infnorm(H)
% [Hinf,fmax] = infnorm(H)	 Find the infinity norm of a z-domain transfer function.

% Get a rough idea of the location of the maximum.
N = 129;
w = linspace(0,2*pi,N);	dw = 2*pi/(N-1);
Hval = evalTF(H,exp(j*w));
[Hinf wi] = max(abs(Hval));

if exist('optimset','file')==2 & exist('fminbnd','file')==2
    % Home in using the "fminbnd" function.
    options = optimset('TolX',1e-8,'TolFun',1e-6);
    wmax = fminbnd('nabsH',w(wi)-dw,w(wi)+dw,options,H);
elseif exist('fmin','file')==2
    % Home in using the "fmin" function.
    options=zeros(1,18);
    options(2)=1e-8;
    options(3)=1e-6;
    wmax = fmin('nabsH',w(wi)-dw,w(wi)+dw,options,H);
else
    fprintf(1,'Hinf: Warning. Optimization toolbox functions not found.\n');
    fprintf(1,' The result returned may not be very accurate.\n');
    wmax = w(wi);
end

Hinf = -nabsH(wmax,H);
fmax = wmax/(2*pi);
