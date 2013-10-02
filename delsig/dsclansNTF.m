function H = dsclansNTF(x,order,rmax,Hz)
% Conversion of clans parameters into a NTF.

% Translate x into H.
% I've changed the relationships between (zeta,wn) and x
% in order to guarantee LHP roots of the s-polynomial.
Hp = zeros(1,length(Hz));
odd = rem(order,2);
if odd
    s = -x(1)^2;
    Hp(1) = rmax*(1+s)./(1-s);
end
for i=1+odd:2:order
    zeta = x(i)^2;
    wn = x(i+1)^2;
    s = roots([1 2*zeta*wn wn^2]);
    Hp(i:i+1) = rmax*(1+s)./(1-s);
end

H = zpk(Hz,Hp,1,1);
