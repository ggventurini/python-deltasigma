function [f1,f2,info]=designHBF(fp,delta,debug)
%function [f1,f2,info]=designHBF(fp=0.2,delta=1e-5,debug=0)
%Design a half-band filter which can be realized without general multipliers.
%The filter is a composition of a prototype and sub- filter.
%Input
% fp	The normalized cutoff frequency of the filter. Due to the
%	symmetry imposed by a HBF, the stopband begins at 0.5-fp.
% delta	The absolute value of the deviation of the frequency response from 
%	the ideal values of 1 in the passband and  0 in the stopband.
%
%Output
% f1,f2	The coefficients of the prototype and sub-filters
%	and their canonical-signed digit (csd) representation.
% info	A vector containing the following data (only set when debug=1):
%	complexity	The number of additions per output sample.
%	n1,n2		The length of the f1 and f2 vectors.
%	sbr		The achieved stob-band attenuation (dB).
%	phi		The scaling factor for the F2 filter.

%Handle the input arguments
parameters = ['fp   ';'delta';'debug'];
defaults = [ 0.2 1e-5 0];
for i=1:length(defaults)
    if i>nargin
       eval([parameters(i,:) '=defaults(i);'])
    elseif eval(['any(isnan(' parameters(i,:) ')) | isempty(' parameters(i,:) ')']) 
       eval([parameters(i,:) '=defaults(i);'])
    end
end

% This is a wrapper function
if exist('firpm','file')
    [f1,f2,info]=designHBF7(fp,delta,debug);
else
    [f1,f2,info]=designHBF6(fp,delta,debug);
end

