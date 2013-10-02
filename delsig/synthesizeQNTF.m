function ntf = synthesizeQNTF(order,OSR,f0,NG,ING,n_im)
%ntf = synthesizeQNTF(order=3,OSR=64,f0=0,NG=-60,ING=-20,n_im=floor(order/3))
%Synthesize a noise transfer function for a quadrature delta-sigma modulator.
%  order  order of the modulator
%  OSR    oversampling ratio
%  f0     center frequency (1->fs)
%  NG     in-band noise gain (dB)
%  ING    image-band noise gain (dB)
%  n_im   number of in-band image zeros
%
%ntf is a zpk object containing the zeros and poles of the NTF. See zpk.m
%
% ALPHA VERSION
% This function uses an experimental ad hoc method that is
% neither optimal nor robust.

% Handle the input arguments
parameters = {'order','OSR','f0','NG','ING','n_im'};
defaults = { 4 64 0 -60 -20 []};
for arg_i=1:length(defaults)
    parameter = char(parameters(arg_i));
    if arg_i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{arg_i};'])
    end
end
if isempty(n_im)
    n_im = floor(order/3);
end
debug_it = 0;

if n_im==0
    % Use synthesizeNTF to get an NTF with the specified NG; ignore ING
    f1 = 0.5/OSR;
    x = 1.5;
    lowest_f = Inf;
    dfdx = NaN;
    for itn = 1:20
	ntf = synthesizeNTF(order,OSR,1,x);
	f = dbv(rmsGain(ntf,0,f1)) - NG;
	if debug_it
	    fprintf(1,'x=%.2f f=%.2f\n', x, f);
	end
	if abs(f) < 0.01
	    break;
	end
	if isnan(dfdx)
	    dx = 0.1*sign(f);
	    dfdx = 0;
	else
	    dfdx = (f-f_old)/dx;
	    dx_old = dx;
	    dx = -f/dfdx;
	    if abs(dx) > max(1,2*abs(dx_old))
		dx = sign(dx)*max(1,2*abs(dx_old));
	    end
	    if x+dx <=1	% Hinf must be at least 1
		dx = dx/2;
	    end
	end
	f_old = f;
	x = x + dx;
    end
    if itn==20
	fprintf(1,'Warning: Iteration limit reached. NTF may be poor\n');
    end
    % Rotate the NTF
    z0 = exp(2i*pi*f0);
    ntf = zpk(z0*ntf.z{1}, z0*ntf.p{1}, ntf.k, 1);
else
    n_in = order-n_im;
    f1 = f0-0.5/OSR;
    f2 = f0+0.5/OSR;
    z0 = exp(2i*pi*f0);
    x = [20 20];	% "R" parameters for cheby2()
    lowest_f = Inf;
    dfdx = [NaN NaN];
    for itn = 1:20
	if debug_it
	    fprintf(1,'x=[%.2f %.2f], ', x(1), x(2) );
	end
	[b1 a1] = cheby2(n_in,x(1),1/OSR,'high');
	[b2 a2] = cheby2(n_im,x(2),1/OSR,'high');
	warning('off');
	ntf0 = zpk( [roots(b1)*z0; roots(b2)*conj(z0)], ...
		    [roots(a1)*z0; roots(a2)*conj(z0)],1,1);
	freq = linspace(-0.5,0.5,200);
	m = evalTF(ntf0,exp(2i*pi*freq));
	NG0 = dbv(rmsGain(ntf0,f1,f2));
	ING0 = dbv(rmsGain(ntf0,-f1,-f2));
	if debug_it
	    clf;
	    subplot(121);
	    plotPZ(ntf0);
	    subplot(122);
	    fprintf(1,'NG=%.1f, ING=%.1f\n', NG0, ING0);
	    plot(freq,dbv(m));
	    figureMagic([-0.5,0.5],0.05,2, [-100 30],10,2)
	    hold on;
	    plot([f1 f2],[1 1]*NG0,'k');
	    text(mean([f1 f2]), NG0, sprintf('NG=%.1fdB',NG0),'vert','bot');
	    plot([-f1 -f2],[1 1]*ING0,'k');
	    text(mean(-[f1 f2]), ING0, sprintf('ING=%.1fdB',ING0),'vert','bot');
	    drawnow;
	end
	f = [NG0 ING0] - [NG ING];
	if abs(f) < 0.01
	    break;
	end
	if norm(f)<lowest_f
	    lowest_f = norm(f);
	    best = ntf0;
	end
	if abs(f(1)) > abs(f(2))
	    i=1;	% adjust x(1)
	else
	    i=2;	% adjust x(2)
	end
	if isnan(dfdx(i))
	    dx = sign(f(i));
	    dfdx(i) = 0;
	    dfdx(3-i) = NaN;
	else
	    dfdx(i) = (f(i)-f_old(i))/dx;
	    dfdx(3-i) = NaN;
	    dx = -f(i)/dfdx(i);
	    xnew = x(i) + dx;
	    if xnew < 0.5*x(i)
		dx = -0.5*x(i);
	    elseif xnew > 2*x(i)
		dx = x(i);
	    end
	end
	f_old = f;
	x(i) = x(i) + dx;
    end
    if itn==20
	fprintf(1,'Warning: Iteration limit reached. NTF may be poor\n');
    end
    ntf = best;
end
