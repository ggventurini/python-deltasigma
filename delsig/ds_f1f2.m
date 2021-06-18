function [f1,f2] = ds_f1f2(OSR,f0,complex_flag);
%[f1,f2] = ds_f1f2(OSR=64,f0=0,complex_flag=0);
parameters = {'OSR' 'f0' 'complex_flag'};
defaults = { 64 0 0 };
for arg_i=1:length(defaults)
    parameter = char(parameters(arg_i));
    if arg_i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{arg_i};'])
    end
end

if complex_flag
    f1 = f0-0.5/OSR;
    f2 = f0+0.5/OSR;
else
    if f0>0.25/OSR
	f1 = f0-0.25/OSR;
	f2 = f0+0.25/OSR;
    else
	f1 = 0;
	f2 = 0.5/OSR;
    end
end
