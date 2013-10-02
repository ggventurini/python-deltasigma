function LCplotTF(H,L0,param)
% LCplotTF(H,L0,param)
%Plot some graphs associated with the LC system.
%
%Arguments
% H	The closed-loop noise transfer function (in z).
% L0	The open-loop tf (in s) from u to v.  
% param	A struct containing the following fields (all normalized):
%   l   The inductances in each tank.
%   c	The capacitances in each tank.
%   gu	The transconductance from the u input to each of the tanks.
%	The final term is the voltage gain to the comparator input.
%   gv	The transconductance from the v output to each of the tanks.
%   gw	The inter-stage transconductances. 
%   gx	The gains from the top of each tank resistor to the comparator.
%   rx	The resistances inserted between the output of each interstage 
%	transconductance and the top of each tank.
%   t	A two-element vector containing the start and end times of the feedback
%	waveform.
%   n	The number of tanks.
%   OSR	The oversampling ratio
%   f0	The center frequency

n = param.n;
OSR = param.OSR;

f0 = param.f0;
fb = 0.5/OSR;

clf;
subplot(321);	% Parameter values
axis([0 20 0 8]);
axis('off')
y = 8;
s=sprintf('g_u='); s=[s sprintf('%g ', param.gu)];
text(0,y,s,'hor','left');	y = y-1;
s=sprintf('g_v='); s=[s sprintf('%g ', param.gv)];
text(0,y,s,'hor','left');	y = y-1;
s=sprintf('g_w='); s=[s sprintf('%g ', param.gw)];
text(0,y,s,'hor','left');	y = y-1;
s=sprintf('g_x='); s=[s sprintf('%g ', param.gx)];
text(0,y,s,'hor','left');	y = y-1;
s=sprintf('r_x='); s=[s sprintf('%g ', param.rx)];
text(0,y,s,'hor','left');	y = y-1;
s = sprintf('t=[%g %g], f_0=fs/%.0f, OSR=%g', param.t(1), param.t(2), ...
 1/f0, param.OSR);
text(0,y,s,'hor','left');	y = y-1;
s=sprintf('l='); s=[s sprintf('%g ', param.l)];
text(0,y,s,'hor','left');	y = y-1;
s=sprintf('c='); s=[s sprintf('%g ', param.c)];
text(0,y,s,'hor','left');	y = y-1;
title('Parameter values');

subplot(322);	% Pole-zero plot
plotPZ(H);
axis('off')
title('Pole-zero plot');
s = sprintf('rmax = %5.3f',max(abs(H.p{1})));
text(0,0,s,'hor','cen','vert','bot');

f = linspace(0,0.5,300);
z = exp(2*pi*j*f);
subplot(323);	% NTF full band
rplot( f, dbv(evalTF(H,z)), [-80 10] );
s = sprintf('Hinf = %4.2f', infnorm(H));
text(0.2,-10,s);
grid on
title('NTF response');
xlabel('normalized frequency');

subplot(324);	% STF full band
resp = dbv(evalTFP(L0,H,f));
rplot( f, resp, [-50 ceil(max(resp)/10)*10] );
pb=find(abs(resp)<3);
fl=f(pb(1)); fu = f(pb(length(pb)));
text(0.25,-5,'3dB frequencies', 'hor', 'cen');
s = sprintf('= %5.3f, %5.3f *fs', fl, fu);
text(0.25,-15,s,'hor','cen');
grid on
title('STF response');
xlabel('1->fs');

fl = f0 - 1.5*fb/2;
fu = f0 + 1.5*fb/2;
f = linspace(fl,fu,100);
z = exp(2*pi*j*f);
subplot(325);	% NTF passband
N=100;
f1 = f0 - fb/2;	f2 = f0 + fb/2;
gain = dbv( rmsGain( H, f1, f2 ) );
ord = ([f1 f2]-f0)/(fb/2);
low = floor(gain/10)*10 -20;
rplot( 2*(f-f0)/fb, dbv(evalTF(H,z)), [low low+40] );
grid on
hold on
plot(ord,[gain gain],'-', ord,[gain gain],'o')
s = sprintf('in-band noise= %3.0fdB', gain-dbp(3*OSR));
text( 0, gain+5, s ,'hor','cen');
xlabel('normalized frequency offset')
hold off

subplot(326);	% STF passband
rplot( (f-f0)/(fb/2), dbv(evalTFP(L0,H,f)), [-3 3]);
grid on
gain1 = dbv( evalTFP(L0,H,(f0-fb/2)) );
gain2 = dbv( evalTFP(L0,H,(f0+fb/2)) );
text( 0, 1.5, 'gain at passband edges =', 'hor', 'center' );
s = sprintf( '%.2fdB, %.2fdB', gain1, gain2 );
text( 0, 0.5, s, 'hor', 'center' );
xlabel('1 \rightarrow f0+fb/2')
return


function [xr,yr]=rplot(x,y,yrange,fmt)
%function rplot(x,y,yrange,fmt)
%plot y vs x, restricted to a certain range.
n = length(x);
if nargin <4
    fmt = '-';
end

%deal with column vectors only
if(size(x,2)~=1)
    x=x';
end
if(size(y,2)~=1)
    y=y';
end

where = (y>yrange(2)) - (y<yrange(1));
above = 1;
inside = 0;
below = -1;

xr=[];	 yr=[];	 i=1;
while i<=n
    tmp = find(where(i:n)~=where(i));
    if isempty(tmp)
	tmp=n+1;
    else
	tmp = tmp(1)+i-1;
    end;
    if where(i)==inside		% copy the data over
	xr = [xr; x(i:tmp-1)];
	yr = [yr; y(i:tmp-1)];
	if( tmp<=n )		% add a point to the data
	    if(where(tmp) == above)
		xr = [xr; x(tmp) + (x(tmp)-x(tmp-1))/(y(tmp)-y(tmp-1))*(yrange(2)-y(tmp))];
		yr = [yr; yrange(2)];
	    else
		xr = [xr; x(tmp) + (x(tmp)-x(tmp-1))/(y(tmp)-y(tmp-1))*(yrange(1)-y(tmp))];
		yr = [yr; yrange(1)];
	    end
	end
    elseif where(i)==above
	if( tmp<=n )		% add a point to the data
	    xr = [xr; x(tmp) + (x(tmp)-x(tmp-1))/(y(tmp)-y(tmp-1))*(yrange(2)-y(tmp))];
	    yr = [yr; yrange(2)];
	    if( where(tmp)==below )	% add another point to the data
		xr = [xr; x(tmp) + (x(tmp)-x(tmp-1))/(y(tmp)-y(tmp-1))*(yrange(1)-y(tmp))];
		yr = [yr; yrange(1)];
	    end
	end
    else
	if( tmp<=n )
	    xr = [xr; x(tmp) + (x(tmp)-x(tmp-1))/(y(tmp)-y(tmp-1))*(yrange(1)-y(tmp))];
	    yr = [yr; yrange(1)];
	    if( where(tmp)==below )	% add another point to the data
		xr = [xr; x(tmp) + (x(tmp)-x(tmp-1))/(y(tmp)-y(tmp-1))*(yrange(2)-y(tmp))];
		yr = [yr; yrange(2)];
	    end
	end
    end
    i = tmp;
end
plot(xr,yr,fmt);
axis([min(x) max(x) yrange(1) yrange(2)]);
return
