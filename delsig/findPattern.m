function [data, snr] = findPattern(N,OSR,ntf,ftest,Atest,f0,nlev,quadrature,dbg)
% [data snr] = findPattern(N=1024,OSR=64,ntf,ftest,Atest,f0,nlev,quadrature,dbg)
% findPattern.m 	Create a length-N bit-stream
% possessing good spectral properties when repeated.
%
% The procedure is:
% Try a bunch of signals. For each signal, calculate the noise power
% induced by resetting the state back to the value N samples earlier.
% This calculation is facilitated by pre-computation of the in-band
% portion of the FFTs of the impulse responses from each integrator's
% input to the output.

% Handle the input arguments
parameters = {'N' 'OSR' 'ntf' 'ftest' 'Atest' 'f0' 'nlev' 'quadrature' 'dbg'};
defaults = { 2^13 64 NaN NaN 0.5 0 2 0 0 };
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin || ( eval(['isnumeric(' parameter ') '])  &&  ...
            eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end
if isnumeric(ntf) && isnan(ntf)
    ntf = synthesizeNTF(6,OSR,1,(nlev+1)/2,f0);
end
if isnan(ftest)
    if quadrature
        if f0==0
            ftest = 0.75/OSR;	% Default ftest is 75% of pb width
        else
            ftest = f0 + 0.75/(2*OSR);
        end
    else
        if f0==0
            ftest = 0.75/(2*OSR);	% Default ftest is 75% of pb width
        else
            ftest = f0 + 0.75/(4*OSR);
        end
    end
end
fbin = round(N*ftest);
ftest = fbin/N;
wtest = 2*pi*ftest;
Nsim = N+1000;			% +1000 is arbitrary

if quadrature
    ABCD = realizeQNTF(ntf,'FB');
    ABCD = mapQtoR(ABCD);
else
    form = 'CRFB';
    [a,g,b,c] = realizeNTF(ntf,form);
    ABCD = stuffABCD(a,g,b,c,form);
end

if quadrature
    if f0==0
        ibb = 0:ceil(N/OSR); % in-band bins
    else
        ibb = floor(N*(f0-1/(2*OSR))):ceil(N*(f0+1/(2*OSR)));
    end
else
    if f0==0
        ibb = 0:ceil(N/(2*OSR)); % in-band bins
    else
        ibb = floor(N*(f0-1/(4*OSR))):ceil(N*(f0+1/(4*OSR)));
    end
end
ibb = setdiff(ibb,fbin);

% Try a bunch of phases and slightly different amplitudes
if length(Atest)<length(wtest)
    Atest(end+1:length(wtest)) = Atest(end);
end
Atest0 = Atest*(nlev-1);
best_fom = Inf;
for itn = 1:50			% 50 is arbitrary
    Atest = Atest0 + 0.01*randn(1,1);
    phi = 2*pi*rand(size(wtest));
    u = zeros(1,Nsim+1);
    for k = 1:length(wtest);
        if quadrature
            u = u + Atest(k) * exp( 1i*( wtest(k)*(0:Nsim) + phi(k) ) );
        else
            u = u + Atest(k) * sin( wtest(k)*(0:Nsim) + phi(k) );
        end
    end
    if quadrature
        % [v xn] = simulateQDSM(u,ABCD,nlev);
        u = [real(u); imag(u)];
        v = simulateDSM(u,ABCD,nlev*[1 1] );
        v = v(1,:) + 1i*v(2,:);
    else
        v = simulateDSM(u,ABCD,nlev);
    end
    for i=0:Nsim-N
        spec = fft(v(i+1:N+i))/(N/2*(nlev-1));
        fom = dbp( sum(abs(spec(ibb+1)).^2) );
        if fom < best_fom
            best_fom = fom;
            best = [Atest*(nlev-1) phi];
            data = v(i+1:N+i);
            if dbg
                fprintf(1,'itn=%d, fom=%g\n', itn, fom);
                make_plot(data,nlev,ibb,fbin);
            end
        end
    end
end
snr = dbv(best(1)) - best_fom;

if dbg
    make_plot(data,nlev,ibb,fbin);
end

return

function make_plot(data,nlev,ibb,fbin)
    N = length(data);
    V = abs(fft(data)) / (N*(nlev-1)/2);
    snr = dbp( sum(V(fbin+1).^2) / sum(V(ibb+1).^2) );
    f = min(ibb):max(ibb);
    plot( f, dbv(V(f+1)), 'b' );
    msg = sprintf(' SNR = %.1f dB', snr);
    text(max(fbin),-20,msg,'hor','left');
    drawnow;
return
