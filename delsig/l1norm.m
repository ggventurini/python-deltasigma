function h1 = l1norm(H)
% h1 = l1norm(H)	Compute the l1-norm of a z-domain transfer function.

if strcmp(H.form,'zp')
    sys = zpk(H.zeros,H.poles,H.k,1);
elseif strcmp(H.form,'coeff')
    sys = tf(H.num,H.den,1);
else
    fprintf(2,'%s: TF form %s is not recognized.\n', mfilename, H.form)
end

h1 = sum(abs(impulse(sys,100)));
