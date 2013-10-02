function f = ds_freq(osr,f0,quadrature)
% f = ds_freq(osr=64,f0=0,quadrature=0)    Frequency vector suitable for plotting the frequency response of an NTF
parameters = {'osr' 'f0' 'quadrature'};
defaults = { 64 0 0 };
for arg_ii=1:length(defaults)
    parameter = parameters{arg_ii};
    if arg_ii>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{arg_ii};'])
    end
end

if quadrature
    f_left = -0.5;
    f_special = [f0 -f0];
else
    f_left = 0;
    f_special = f0;
end
f = linspace(f_left,0.5,100);
% Use finer spacing in the vicinity of the passband
for fx = f_special
    f1 = max(f_left, fx-1/osr);
    f2 = min(0.5,fx+2/osr);
    f( f<=f2 & f>=f1 ) = [];
    f = sort( [f linspace(f1,f2)] );
end
