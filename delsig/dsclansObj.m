function [f,g]=dsclansObj(x,order,OSR,Q,rmax,Hz)
% Objective function for clans.m 
% f is the mean-square value of H over the passband
% g = (||h||_1 - Q)^2

% Translate x into H.
H = dsclansNTF(x,order,rmax,Hz);

% Compute f and g
z = exp(1i*linspace(0,pi/OSR));
f = norm(evalTF(H,z))^2;
g = (sum(abs(impulse(H,100))) -1 - Q)^2;
