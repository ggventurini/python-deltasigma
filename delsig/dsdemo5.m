% Demonstrate the simulateESL function
echo off;
if exist('LiveDemo','var') == 0
    LiveDemo=0;
end
clc;
fprintf(1,'\t\tMismatch-Shaping Unit-Element DAC\n\n');

% Specify the modulator NTF, the mismatch-shaping TF, and the number of elements
echo on;
ntf = synthesizeNTF(3,[],[],4);
M = 16;
sigma_d = 0.01;		% 1% mismatch
mtf1 = zpk(1,0,1,1);	%First-order shaping
echo off;	
A = 1/sqrt(2);		% Test tone amplitude, relative to full-scale.
f = 0.33;			% Test tone frequency, relative to fB. 
					% (Will be adjusted to be an fft bin)

if LiveDemo
    N = 2^12;
else
    N = 2^14;
end
fin = round(f*N/(2*12));
w = (2*pi/N)*fin;
echo on;
u = M*A*sin(w*[0:N-1]);
v = simulateDSM(u,ntf,M+1);	% M unit elements requires an M+1-level quant.
v = (v+M)/2;				% scale v to [0,M]
sv0 = thermometer(v,M);
sv1 = simulateESL(v,mtf1,M);
echo off

figure(1); clf
T = 20;
subplot(211);
plotUsage(sv0(:,1:T));
set(gcf,'NumberTitle','off'); 
set(gcf,'Name','Element Usage');
title('Thermometer-Coding')
subplot(212);
plotUsage(sv1(:,1:T));
title('First-Order Shaping');
if LiveDemo
    set(1,'position',[9 204 330 525]);
    changeFig(18,.5,1);
    pause
end

ideal = v;

% DAC element values
e_d = randn(M,1); 
e_d = e_d - mean(e_d);
e_d = sigma_d * e_d/std(e_d);
ue = 1 + e_d;

% Convert v to analog form, assuming no shaping
dv0 = ue' * sv0;

% Convert sv to analog form 
dv1 = ue' * sv1;

window = ds_hann(N);
spec = fft(ideal.*window)/(M*N/8);
spec0 = fft(dv0.*window)/(M*N/8);
spec1 = fft(dv1.*window)/(M*N/8);

figure(2); clf
set(gcf,'NumberTitle','off'); 
set(gcf,'Name','Spectra');
plotSpectrum(spec0,fin,'r');
hold on;
plotSpectrum(spec1,fin,'b');
plotSpectrum(spec,fin,'g');
axis([1e-3 0.5 -200 0]);
x1 = 2e-3; x2=1e-2; y0=-180; dy=dbv(x2/x1); y3=y0+3*dy;
plot([x1 x2 x2 x1],[y0 y0 y3 y0],'k')
text(x2, (y0+y3)/2,' 60 dB/decade')
hold off;
grid;
ylabel('PSD');
xlabel('Normalized Frequency');
if LiveDemo
    figure(1);
    set(1,'position',[9 427 200 300]);
    changeFig;
    subplot(211)
    xlabel('')
    figure(2);
    set(2,'position',[237 288 553 439]);
    changeFig(18,2,12);
    legend('thermom.','1^{st}-Order','Ideal');
    set(gcf,'NumberTitle','off'); 
    set(gcf,'Name','Output Spectra');
	pause
    changeFig;
    legend;
    set(2,'position',[238 427 325 300]);
end
legend('thermometer','rotation','ideal DAC');
fprintf(1,'Paused.\n');
pause

%Now repeat the above for second-order shaping
echo on;
mtf2 = zpk([ 1 1 ], [ 0 0 ], 1, 1);	%Second-order shaping
sv2 = simulateESL(v,mtf2,M);
echo off;

figure(1); clf
T = 20;
subplot(211);
plotUsage(sv1(:,1:T));
title('First-Order Shaping');
subplot(212);
plotUsage(sv2(:,1:T));
title('Second-Order Shaping');
if LiveDemo
    set(1,'position',[9 204 330 525]);
    changeFig(18,.5,1);
    pause
end

dv2 = ue' * sv2;
spec2 = fft(dv2.*window)/(M*N/8);

figure(2); clf
plotSpectrum(spec1,fin,'r');
hold on;
plotSpectrum(spec2,fin,'b');
plotSpectrum(spec,fin,'g');
axis([1e-3 0.5 -200 0]);
x1 = 2e-3; x2=1e-2; y0=-180; dy=dbv(x2/x1); y3=y0+3*dy;
plot([x1 x2 x2 x1],[y0 y0 y3 y0],'k')
text(x2, (y0+y3)/2,' 60 dB/decade')
legend('1^{st}-Order','2^{nd}-Order','Ideal');
hold off;
grid;
xlabel('Normalized Frequency');
ylabel('PSD');
if LiveDemo
    figure(1);
    set(1,'position',[9 427 200 300]);
    changeFig;
    subplot(211)
    xlabel('')
    figure(2);
    set(2,'position',[237 288 553 439]);
    changeFig(18,2,12);
    legend('1^{st}-Order','2^{nd}-Order','Ideal');
	pause
    changeFig;
    legend('1^{st}-Order','2^{nd}-Order','Ideal');
    set(2,'position',[238 427 325 300]);
end

if 0 
    % Plot everything
    figure(1); clf
    T = 20;
    subplot(311);
    plotUsage(thermometer(v(1:T),M));
    set(gcf,'NumberTitle','off'); 
    set(gcf,'Name','Element Usage');
    title('Thermometer-Coding')
    subplot(312);
    plotUsage(sv1(:,1:T));
    title('First-Order Shaping');
    subplot(313);
    plotUsage(sv2(:,1:T));
    title('Second-Order Shaping');
    %printmif(fullfile('MIF','elementUsage'), [3 6], 'Helvetica10')
    figure(2); clf
    plotSpectrum(spec0,fin,'r');
    hold on;
    plotSpectrum(spec1,fin,'r');
    plotSpectrum(spec2,fin,'b');
    plotSpectrum(spec,fin,'g');
    axis([1e-3 0.5 -200 0]);
    x1 = 2e-3; x2=1e-2; y0=-180; dy=dbv(x2/x1); y3=y0+3*dy;
    plot([x1 x2 x2 x1],[y0 y0 y3 y0],'k')
    text(x2, (y0+y3)/2,' 60 dB/decade')
    legend('Thermometer','1^{st}-Order','2^{nd}-Order','Ideal');
    hold off;
    grid;
    xlabel('Normalized Frequency');
    ylabel('PSD');
    printmif(fullfile('MIF','mmSpetcra'), [3 3], 'Helvetica10')
end

if 0 % Quadrature example
    order = 4;
    osr = 32;
    M = 8;
    f0 = 1/16;
    quadrature = 1;
    [f1 f2] = ds_f1f2(osr,f0,quadrature);
    ntf = synthesizeQNTF(order,osr,f0,-50,-10);
    Amp = undbv(-3);  % Test tone amplitude, relative to full-scale.
    f = 0.3;          % Test tone frequency offset from f0, relative to bw. 
                      % (Will be adjusted to be an fft bin)
    N = 2^12;
    f1_bin = round(f1*N);
    f2_bin = round(f2*N);
    fin = round(((1-f)/2*f1 + (f+1)/2*f2)*N);

    u = Amp*M*exp((2i*pi/N)*fin*[0:N-1]);
    v = simulateQDSM(u,ntf,M+1);
    sv0 = [2i*thermometer((imag(v)+M)/2,M)-1i; 
           2*thermometer((real(v)+M)/2,M)-1];
    mtf1 = zpk(exp(2i*pi*f0),0,1,1);	%First-order complex shaping
    sv1 = simulateQESL(v,mtf1,M);
    figure(1); clf
    subplot(121);
    plotUsage(sv0(:,1:20));
    subplot(122);
    plotUsage(sv1(:,1:20));
    %printmif(fullfile('MIF','QUsage'), [32 4], 'Helvetica8')    
    ue = 1 + 0.01*randn(2*M,1);
    dv0 = ue' * sv0;
    dv1 = ue' * sv1;
    window = ds_hann(N);
    Nsm = 32;
    NBW = 1.5/N * Nsm/2 / 1.5;  % Approximate correction for smoothing; amplitude correction is also approximate
    spec =  fft(v.*window)/(M*N/2);
    spec0 = fft(dv0.*window)/(M*N/2);
    spec1 = fft(dv1.*window)/(M*N/2);
    spec_sm = circ_smooth(abs(fftshift(spec)).^2,Nsm)* Nsm/3;
    spec0_sm = circ_smooth(abs(fftshift(spec0)).^2,Nsm)* Nsm/3;
    spec1_sm = circ_smooth(abs(fftshift(spec1)).^2,Nsm)* Nsm/3;
    freq = linspace(-0.5,0.5,N+1); freq(end)=[];
    figure(2); clf
	plot(freq,dbp(spec_sm),'k','Linewidth',1);
    hold on;
	plot(freq,dbp(spec0_sm),'r','Linewidth',1);
	plot(freq,dbp(spec1_sm),'b','Linewidth',1);
    snr = calculateSNR(spec(f1_bin+1:f2_bin+1),fin-f1_bin);
    snr0 = calculateSNR(spec0(f1_bin+1:f2_bin+1),fin-f1_bin);
    snr1 = calculateSNR(spec1(f1_bin+1:f2_bin+1),fin-f1_bin);
    msg = sprintf('Ideal: SNR = %.0fdB', snr);
    msg0 = sprintf('No MS: SNR = %.0fdB', snr0);
    msg1 = sprintf('1^{st}-order MS: SNR = %.0fdB', snr1);
    text(0.5,-95,sprintf('NBW=%.0e ',NBW),'Hor','right');
    figureMagic([-0.5 0.5],1/16,4, [-100 0],10,2,[6 3],'Spectra');
    xlabel('frequency')
    ylabel('PSD (dBFS/NBW)');
    legend(msg,msg0,msg1);
    printmif(fullfile('MIF','QMS'), [4 2], 'Helvetica8')
end
