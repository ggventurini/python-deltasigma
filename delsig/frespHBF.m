function [mag,pbr,sbr] = frespHBF(f,f1,f2,phi,fp,msg)
%[mag pbr sbr] = frespHBF(f,f1,f2,phi,fp,msg)	
%Compute the frequency response, the passband ripple and the stopband ripple 
% of a Saramaki HBF. If msg is non-null, a plot is made.
% fp is the passband edge.
% phi is used by designHBF.

% Handle the input arguments
parameters = {'f' 'f1' 'f2' 'phi' 'fp' 'msg'};
defaults = { NaN NaN NaN 1 0.2 '' };
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end
if isnan(f)
    f = linspace(0,.5,1024);
end
if isstruct(f1)	% Presume that f1 is a {.val,.csd} struct
    f1 = [f1.val];
end
if isstruct(f2)	% Presume that f2 is a {.val,.csd} struct
    f2 = [f2.val];
end
f1 = f1(:); f2 = f2(:);

Npts = length(f);
w = 2*pi*f;
z = exp(j*w);
cos_w = real(z);

n2 = length(f2);
F2 = zeros(size(w));
for i = 1:n2
    F2 = F2 + f2(i)*cos(w*(2*i-1));
end
F2 = F2*2;
mag = evalF1(f1,F2); 

if nargout > 1 | ~isempty(msg)
    passband = 1:floor(2*fp*(Npts-1) +1);
    stopband = Npts + 1 - passband;
    pbr = max( abs( abs(mag(passband)) -1 ) );
    sbr = max( abs(mag(stopband)) );
end

if ~isempty(msg)
    clf
    subplot(211)
    F1 = evalF0(f1,z,phi); 
    hndl(1) = plot(f,abs(F1),'--'); hold on;
    hndl(2) = plot(f,phi*abs(F2),':');
    hndl(3) = plot(f, abs(mag),'-'); hold off;
    legend(hndl,'F1','F2','HBF')
    title(msg);
    hold on; grid on
    axis([0 0.5 0 1.1])

    subplot(212)
    plot(f,dbv(mag))
    axis([0 0.5 -150 10])
    grid on

    msg = sprintf( ' pbr=%.1e', pbr );
    text(0.0, -10, msg, 'VerticalAlignment', 'top');
    msg = sprintf( 'sbr=%.0fdB ', dbv(sbr) );
    text(0.5, dbv(sbr), msg, 'HorizontalAlignment', 'right', ...
      'VerticalAlignment', 'bottom');
    drawnow;
    set(gcf,'MenuBar','none');
    set(gcf,'NumberTitle','off');
    set(gcf,'Name','HBF Frequency Response');
end
