function NTF=clans5(order,OSR,Q,rmax,opt)
%Version of clans for MATLAB 5 or lower

% Create the initial guess
NTF = synthesizeNTF(order,OSR,opt,1+Q,0);
Hz = NTF.z{1};
x = zeros(1,order);
odd = rem(order,2);
if odd
    z = NTF.p{1}(1)/rmax;
    if any(abs(z))>1 %project poles outside rmax onto the circle
	z = z./abs(z);
    end
    s = (z-1)./(z+1);
    x(1)= sqrt(-s);
end
for i=odd+1:2:order
    z = NTF.p{1}(i:i+1)/rmax;
    if any(abs(z))>1 %project poles outside rmax onto the circle
	z = z./abs(z);
    end
    s = (z-1)./(z+1);
    coeffs=poly(s);
    wn = sqrt(coeffs(3));
    zeta = coeffs(2)/(2*wn);
    x(i) = sqrt(zeta);
    x(i+1) = sqrt(wn);
end

% Run the optimizer
x=constr('dsclansObj',x,[],[],[],[],order,OSR,Q,rmax,Hz);
NTF = dsclansNTF(x,order,rmax,Hz);
return

