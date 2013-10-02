function ntf = synthesizeNTF(order,osr,opt,H_inf,f0)
%ntf = synthesizeNTF(order=3,osr=64,opt=0,H_inf=1.5,f0=0)
%Synthesize a noise transfer function for a delta-sigma modulator.
%	order =	order of the modulator
%	osr =	oversampling ratio
%	opt =	flag for optimized zeros
%		0 -> not optimized,
%		1 -> optimized, 
%		2 -> optimized with at least one zero at band-center
%		3 -> optimized zeros (Requires MATLAB6 and Optimization Toolbox)
%       [z] -> zero locations in complex form
%	H_inf =	maximum NTF gain
%	f0 =	center frequency (1->fs)
%
%ntf is a zpk object containing the zeros and poles of the NTF. See zpk.m
%
% See also 
%  clans()   "Closed-loop analysis of noise-shaper." An alternative
%            method for selecting NTFs based on the 1-norm of the 
%            impulse response of the NTF
%
%  synthesizeChebyshevNTF()    Select a type-2 highpass Chebyshev NTF.
%            This function does a better job than synthesizeNTF if osr
%            or H_inf is low.

% This is actually a wrapper function which calls either the 
% appropriate version of synthesizeNTF based on the availability
% of the 'fmincon' function from the Optimization Toolbox

% Handle the input arguments
parameters = {'order' 'osr' 'opt' 'H_inf' 'f0'};
defaults = { 3 64 0 1.5 0 };
for arg_i=1:length(defaults)
    parameter = char(parameters(arg_i));
    if arg_i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{arg_i};'])
    end
end
if f0 > 0.5 
    fprintf(1,'Error. f0 must be less than 0.5.\n');
    return;
end
if f0 ~= 0 & f0 < 0.25/osr
    warning('(%s) Creating a lowpass ntf.', mfilename);
    f0 = 0;
end
if f0 ~= 0 & rem(order,2) ~= 0
    fprintf(1,'Error. order must be even for a bandpass modulator.\n');
    return;
end

if length(opt)>1 & length(opt)~=order
    fprintf(1,'The opt vector must be of length %d(=order).\n', order);
    return;
end


if exist('fmincon','file')
    ntf = synthesizeNTF1(order,osr,opt,H_inf,f0);
else
    ntf = synthesizeNTF0(order,osr,opt,H_inf,f0);
end
