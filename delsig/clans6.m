function NTF=clans6(order,OSR,Q,rmax,opt)
%Version of clans for MATLAB >=6

% Create the initial guess
NTF = synthesizeNTF(order,OSR,opt,1+Q,0);
Hz = NTF.z{1};
x = zeros(1,order);
odd = rem(order,2);
poles = cplxpair(NTF.p{1}); poles=poles(end:-1:1);
if odd
    z = poles(1)/rmax;
    if any(abs(z))>1 %project poles outside rmax onto the circle
	z = z./abs(z);
    end
    s = (z-1)./(z+1);
    x(1)= sqrt(-s);
end
for i=odd+1:2:order
    z = poles(i:i+1)/rmax;
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
options = optimset('TolX',1e-6, 'TolFun',1e-6, 'TolCon',1e-6, 'MaxIter',1000 );
options = optimset(options,'Display','off');
options = optimset(options,'Diagnostics','off');
options = optimset(options,'LargeScale','off');
x = fmincon(@(x)dsclansObj6a(x,order,OSR,Q,rmax,Hz),x, ...
[],[],[],[],[],[],@(x)dsclansObj6b(x,order,OSR,Q,rmax,Hz), options );

NTF = dsclansNTF(x,order,rmax,Hz);
return

function f=dsclansObj6a(x,order,OSR,Q,rmax,Hz)
% Objective function for clans; Optimization Toolbox version >= 6
% f is the magnitude of H at the band-edge
H = dsclansNTF(x,order,rmax,Hz);
f = abs(evalTF(H,exp(1i*pi/OSR)));
return

function [g,g_eq]=dsclansObj6b(x,order,OSR,Q,rmax,Hz)
% Constraint function for clans; Optimization Toolbox version >= 6
% g =||h||_1 - Q
H = dsclansNTF(x,order,rmax,Hz);
g = sum(abs(impulse(H,100))) -1 - Q;
g_eq = [];
return
