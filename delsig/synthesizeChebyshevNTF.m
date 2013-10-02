function ntf = synthesizeChebyshevNTF(order,OSR,opt,H_inf,f0)
%ntf = synthesizeChebyshevNTF(order=3,OSR=64,opt,H_inf=1.5,f0=0)
%Synthesize a noise transfer function for a delta-sigma modulator.
%The NTF is a type-2 highpass Chebyshev function
%	order =	order of the modulator
%	OSR =	oversampling ratio
%	H_inf =	maximum NTF gain
%	f0 =	center frequency (1->fs)
%
%ntf is a zpk object containing the zeros and poles of the NTF. See zpk.m
%

% Handle the input arguments
parameters = {'order' 'OSR' 'opt' 'H_inf' 'f0'};
defaults = { 3 64 0 1.5 0 };
for arg_i=1:length(defaults)
    parameter = char(parameters(arg_i));
    if arg_i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{arg_i};'])
    end
end
if( f0 ~= 0 )
    if( rem(order,2) ~= 0)
	fprintf(1,'order must be even for a bandpass modulator.\n');
	return;
    else
	[f1 f2] = ds_f1f2(OSR,f0);
	f1f2 = [f1 f2];
    end
end

% Iteratively solve for the attenuation spec (x) which yields the desired H_inf
x_min = 0;
x_max = 300;
dx_max = 10;
ftol = 1e-6;
xtol = 1e-6;
x = 60;		% Initial guess
itn_limit = 10;
converged = 0;
for itn = 1:itn_limit
    if( f0 == 0 )
	[z p k] = cheby2(order, x, 1/OSR, 'high');
    else
	[z p k] = cheby2(order/2, x, 2*f1f2, 'stop');
    end
    f = 1/k - H_inf;
    % fprintf(1,'%9.6f %9.6f\n', x, f);
    if f>0	% x is too big
	x_max = x;
    else	% x is too small
	x_min = x;
    end
    if itn==1
	dx = -dx_max * sign(f);
    else
	df = f - f_p;
	if abs(df) < ftol 
	    converged = 1;
	    break;
	end
	dx = -f * dx/df;
    end
    if converged
	break;
    end
    x_p = x;
    f_p = f;
    x = max( [x_min min([x+dx,x_max])] );
    dx = x - x_p;
    if( abs(dx) < xtol )
	break;
    end
end
ntf = zpk(z,p,1,1);
