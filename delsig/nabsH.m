function y=nabsH(w,H)
% nabsH computes the negative of the absolute value of H 
%       at the specified frequency on the unit circle. 
%
%       This function is used by infnorm.m.
z = exp(i*w);
y = -abs(evalTF(H,z));
