function [f,g]=dsclansObj(x,order,OSR,Q,rmax,Hz)
% Objective function for clans.m 
% f is the magnitude of H at the band-edge
% g =||h||_1 - Q

% Translate x into H.
H = dsclansNTF(x,order,rmax,Hz);

% Compute f and g
I = sqrt(-1);
f = abs(evalTF(H,exp(I*pi/OSR)));
g = sum(abs(impulse(H,100))) -1 - Q;
