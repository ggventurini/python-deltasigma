function [snr,amp,k0,k1,sigma_e2] = predictSNR(ntf,R,amp,f0)
%[snr,amp,k0,k1,sigma_e2] = predictSNR(ntf,R=64,amp=...,f0=0)
%Predict the SNR curve of a binary delta-sigma modulator by using
%the describing function method of Ardalan and Paulos.
%
%The modulator is specified by a noise transfer function (ntf).
%The band of interest is defined by the oversampling ratio (R)
%and the center frequency (f0).
%The input signal is characterized by the amp vector.
%amp defaults to [-120 -110...-20 -15 -10 -9 -8 ... 0]dB, where 0 dB means
%a full-scale (peak value = 1) sine wave. 
%
%The algorithm assumes that the amp vector is sorted in increasing order;
%once instability is detected, the remaining SNR values are set to -Inf.
%
% Output:
%snr 	a vector of SNR values (in dB)
%amp	a vector of amplitudes (in dB)
%k0	the quantizer signal gain
%k1	the quantizer noise gain
%sigma_e2 the power of the quantizer noise (not in dB)
%
%The describing function method of A&P assumes that the quantizer processes
%signal and noise components separately. The quantizer is modelled as two
%(not necessarily equal) linear gains, k0 and k1, and an additive white
%gaussian noise source of power sigma_e2. k0, k1 and sigma_e2 are calculated
%as functions of the input.
%The modulator's loop filter is assumed to have nearly infinite gain at
%the test frequency.
%
% Future versions may accommodate STFs.

% Handle the input arguments
if nargin<1
    error('Insufficient arguments');
end
parameters = {'ntf';'R';'amp';'f0'};
defaults = {NaN, 64, [-120:10:-20 -15 -10:0], 0};
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end

Nb = 100;
if f0==0
    band_of_interest = linspace(0,pi/R,Nb);
else
    band_of_interest = linspace(2*pi*(f0-0.25/R),2*pi*(f0+0.25/R),Nb);
    XTAB = linspace(-2,0,21);
    % The following code was used to create the YTAB matrix
%    YTAB = [];
%    for xi=XTAB
%	YTAB = [YTAB; hypergeo(0.5,1,xi) hypergeo(0.5,2,xi)];
%    end
    YTAB = [0.46575960516930   0.67366999387741
	    0.47904652357101   0.68426650762558
	    0.49316295981407   0.69527947902679
	    0.50817364454269   0.70673173666000
	    0.52414894104004   0.71864765882492
	    0.54116523265839   0.73105299472809
	    0.55930554866791   0.74397552013397
	    0.57866013050079   0.75744456052780
	    0.59932720661163   0.77149158716202
	    0.62141352891922   0.78615015745163
	    0.64503526687622   0.80145609378815
	    0.67031890153885   0.81744754314423
	    0.69740217924118   0.83416539430618
	    0.72643494606018   0.85165339708328
	    0.75758063793182   0.86995816230774
	    0.79101717472076   0.88912981748581
	    0.82693856954575   0.90922164916992
	    0.86555624008179   0.93029111623764
	    0.90710091590881   0.95239937305450
	    0.95182400941849   0.97561222314835
	    1.00000000000000   1.00000000000000];
end

[num,den] = zp2tf(ntf.z{:},ntf.p{:},1);
num1 = num - den;

N = length(amp);
snr = zeros(1,N)-Inf;
k0 = zeros(1,N);
k1 = zeros(1,N);
sigma_e2 = zeros(1,N);

u = 10.^(amp/20);

Nimp = 100;
unstable = 0;
for n=1:N;
    % Calculate sigma_e2
    if f0==0
	erfinvu = erfinv(u(n));
	sigma_e2(n) = 1 - u(n)^2 - 2/pi*exp(-2*erfinvu^2);
    else % Sinusoidal input.
	% Solve sqrt(pi)*u/2 = rho * hypergeo(0.5,2,-rho^2);
	% Formulate as solve f(rho)=0, f = rho*M(0.5,2,-rho^2)-K
	% and use the secant method.
	K = 0.5 * sqrt(pi) * u(n);
	if n==1
	    rho = u(n)^2;	% Initial guess; otherwise use previous value.
	    fprime = 1;
	end
	drho = 1;
	for itn = 1:20
	    m0 = interp1(XTAB,YTAB(:,2),-rho^2,'*cubic');
	    f = rho*m0 - K;
	    if( itn >1 )
		fprime = max((f-f_prev)/drho,0.5);	%Secant approx.
	    end
	    if abs(f) < 1e-8; break; end	%!Converged
	    drho = -f/fprime;
	    if abs(drho) > 0.2; drho = sign(drho)*0.2; end
	    if abs(drho) < 1e-6; break; end	%!Converged
	    rho = rho + drho;
	    f_prev = f;
	end
	m1 = interp1(XTAB,YTAB(:,1),-rho^2,'*cubic');
	sigma_e2(n) = 1 - u(n)^2/2 - 2/pi*m1^2;
    end

    % Iterate to solve for k1 and sigma_1.
    % Using one of MATLAB's nonlinear equation solvers would be more efficient,
    % but this function code would then require the optimization toolbox.
    % !Future work: put in 2-D BFGS code.
    if n>1
	k1(n) = k1(n-1); % Use the previous value of k1 as the initial guess.
    else
	k1(n) = 1.2;
    end
    k1_prev = 0;
    itn = 0;
    if f0==0
	k1sigma1 = sqrt(2/pi)*exp(-erfinvu^2);
    else
	k1sigma1 = sqrt(2/pi)*m1;
    end
    while abs(k1(n)-k1_prev) > 1e-6*(1+k1(n)) & itn < 100
	% Create the function: H_hat = L1/(1-k1*L1)=(H-1)/(H*(1-k1)+k1).
	den1 = (1-k1(n))*num + den*k1(n);
	% Calculate pGain, the square of the 2-norm of H_hat.
	[pGain Nimp] = powerGain(num1,den1,Nimp);
	if isinf(pGain)
	    unstable = 1;
	    break
	end

	sigma_1 = sqrt(pGain*sigma_e2(n));
	k1_prev = k1(n);
	k1(n) = k1sigma1/sigma_1;
	itn = itn+1;
    end
    if unstable
	break
    end

    if f0==0 
	y0 = sqrt(2)*erfinvu*sigma_1;
	k0(n) = u(n)/y0;
    else
	k0(n) = sqrt(2/pi)*m0/sigma_1;
    end

    h = freqz(num, (1-k1(n))*num + k1(n)*den, band_of_interest);
    % For both DC and sine wave inputs, use u^2/2 as the signal 
    % power since true DC measurements are usually impossible.
    snr(n) = dbp( 0.5*u(n)^2 / (sum(h.^2)/(R*Nb)*sigma_e2(n)) );
end

function [pGain, Nimp] = powerGain(num,den,Nimp0)
%[pGain Nimp] = powerGain(num,den,Nimp0=100) Calculate the power gain
%of a TF given in coefficient form.
%Nimp is the recommended number of impulse response samples for use
%in future calls and Nimp0 is the suggested number to use.
if nargin<3
    Nimp0=100;
end
Nimp = Nimp0;
unstable = 0;

sys = tf(num,den,1);
imp = impulse(sys,Nimp);
if( sum(abs(imp(Nimp-10:Nimp))) < 1e-8 & Nimp > 50) % Use fewer samples.
    Nimp = round(Nimp/1.3);
else
    while( sum(abs(imp(Nimp-10:Nimp))) > 1e-6 )
	Nimp = Nimp*2;
	imp = impulse(sys,Nimp);
	if sum(abs(imp(Nimp-10:Nimp))) >= 50 | Nimp >= 1e4
	    % H is close to being unstable
	    unstable = 1;
	    break;
	end
    end
end
if unstable==0
    pGain = sum(imp.^2);
else
    pGain = Inf;
end
