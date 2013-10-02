function NTF=clans(order,OSR,Q,rmax,opt)
%NTF = clans(order=4,OSR=64,Q=5,rmax=0.95,opt=0)	Optimal NTF design
%for a multi-bit modulator.
%CLANS = "closed-loop analysis of noise-shapers,"
%and was originally developed by J.G. Kenney and L.R. Carley.

% Handle the input arguments
parameters = {'order';'OSR';'Q';'rmax';'opt'};
defaults = [ 4 64 5 0.95 0 ];
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults(i);'])
    end
end

if exist('fmincon','file')==2		% Optimization Toolbox >=6
    NTF=clans6(order,OSR,Q,rmax,opt);
elseif exist('constr','file')==2	% Optimization Toolbox version < 6
    NTF=clans5(order,OSR,Q,rmax,opt);
else
	error('CLANS needs the optimization toolbox.');
end

