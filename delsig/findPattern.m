function [data, snr] = findPattern(N,OSR,NTF,ftest,Atest,f0,nlev,quadrature,debug)
% [data snr] = findPattern(N=1024,OSR=64,NTF,ftest,Atest,f0,nlev,quadrature,debug)
% findPattern.m 	Create a length-N bit-stream 
% possessing good spectral properties when repeated.
%
% Iniitially created March 16 2000 for Jon Baldwin of ADI/CTS
% March 17 2000: Converted to function form, but not tested
% March 27 2001: Corrected some errors, generalized to bandpass, added two-tone support
% June 13 2003: Generalized to quadrature
%
%
% The procedure is:
% Try a bunch of signals. For each signal, calculate the noise power
% induced by resetting the state back to the value N samples earlier.
% This calculation is facilitated by pre-computation of the in-band
% portion of the FFTs of the impulse responses from each integrator's
% input to the output.

% Handle the input arguments
parameters = {'N' 'OSR' 'NTF' 'ftest' 'Atest' 'f0' 'nlev' 'quadrature' 'debug'};
defaults = { 2^13 64 NaN NaN 0.5 0 2 0 0 };
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end
if isnumeric(NTF) & isnan(NTF)
    NTF = synthesizeNTF(6,OSR,1,(nlev+1)/2,f0);
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
Nsim = N + 2000;			% +1000 is arbitrary

% Calculate the in-band portion of the FFTs of the impulse responses
if quadrature
    ABCD = realizeQNTF(NTF,'FB');
    [A,B,C,D] = partitionABCD(ABCD);
    B2 = B(:,2);
    order = size(A,1);
    D = zeros(1,order);
    A = A + B2*C;
    B = eye(size(A));
    sys = ss(mapQtoR(A), mapQtoR(B), mapQtoR(C), mapQtoR(D), 1);
    imp = impulse(sys,N-1);
    imp = squeeze( imp(:,1,1:2:2*order) +1i*imp(:,2,1:2:2*order) );
    ABCD = mapQtoR(ABCD);
else
    form = 'CRFB';
    [a,g,b,c] = realizeNTF(NTF,form);
    ABCD = stuffABCD(a,g,b,c,form);
    [A,B,C,D] = partitionABCD(ABCD);
    B2 = B(:,2);
    order = size(A,1);
    D = zeros(1,order);
    A = A + B2*C;
    B = eye(size(A));
    sys = ss(A,B,C,D,1);
    imp = impulse(sys,N-1);
    imp = squeeze(imp);
end

if quadrature
    if f0==0
	ibb = [1:ceil(N/OSR+1)]; % in-band bins
    else
	ibb = [floor(N*(f0-1/(2*OSR))):ceil(N*(f0+1/(2*OSR)))] +1;
    end
else
    if f0==0
	ibb = [1:ceil(N/(2*OSR)+1)]; % in-band bins
    else
	ibb = [max(1,floor(N*(f0-1/(4*OSR)))):ceil(N*(f0+1/(4*OSR)))] +1;
    end
end

in_band = zeros(length(ibb),order);
for i = 1:order
    spec = fft(imp(:,i))/(N/2);
    in_band(:,i) = spec(ibb);
end

% Try a bunch of phases and slightly different amplitudes
Atest0 = Atest;
best_fom = Inf;
for i = 1:100			% 50 is arbitrary
    Atest = Atest0 + 0.01*randn(1,1);
    phi = 2*pi*rand(size(wtest));
    u = zeros(1,Nsim+1);
    for k = 1:length(wtest);
	if quadrature
	    u = u +  exp( 1i*( wtest(k)*[0:Nsim] + phi(k) ) );
	else
	    u = u +  sin( wtest(k)*[0:Nsim] + phi(k) );
	end
    end
    u = u*Atest*(nlev-1)/length(wtest);
    if quadrature
	% [v xn] = simulateQDSM(u,ABCD,nlev);
	u = [real(u); imag(u)];
	[v xn] = simulateDSM(u,ABCD,nlev*[1 1] );
	v = v(1,:) + 1i*v(2,:);
	xn = xn(1:2:end,:) + 1i*xn(2:2:end,:);
    else
	[v xn] = simulateDSM(u,ABCD,nlev);
    end
    delta = xn(:,1:Nsim-N) - xn(:,N+1:Nsim);
    spec = in_band * delta;
    [fom start] = min( dbp( sum(abs(spec).^2) ) );
    if fom < best_fom
	best_fom = fom;
	best = [Atest*(nlev-1) phi start];
	data = v(start+1:start+N);
	if debug
	    fprintf(1,'i=%d, fom=%g\n', i, fom);
	end
    end
end
snr = dbv(best(1)) - best_fom;

if 1 | debug
    V0 = fft(data.*ds_hann(N))/(N*(nlev-1)/4);
    snr0 = calculateSNR(V0(ibb),fbin-ibb(1)+1);
    V1 = fft(data)/(N*(nlev-1)/2);
    snr1 = -dbp( sum(abs(V1(ibb).^2)) / sum(abs(V1(fbin+1)).^2) -1 );
    f = ibb-1;
    plot( f,dbv(V0(ibb)), f,dbv(V1(ibb)) )
    msg = sprintf('SNR = %.1fdB (%.1fdB)', snr1, snr0);
    text(f(round(end/2)),-20,msg,'hor','center');
end

return
