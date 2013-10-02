function adc = dsexample1
% Design example for a lowpass modulator.

% Design parameters
order = 5;
osr = 32;
nlev = 2;
f0 = 0;
Hinf = 1.5;
form = 'CRFB';

% Derived parameters
M = nlev-1;

clc
fprintf(1,'\t\t\t%dth-Order Lowpass Example\n\n',order);
fig_pos = { [  10 595  600 215],
            [ 640 595  600 215],
            [  10 345  250 205],
            };

% NTF synthesis and realization
fprintf(1,'Doing NTF synthesis and realization... ');
design_step = 1;
ntf = synthesizeNTF(order,osr,2,Hinf,f0);		% Optimized zero placement
[a,g,b,c] = realizeNTF(ntf,form);
b = [b(1) zeros(1,length(b)-1)];	% Use a single feed-in for the input
ABCD = stuffABCD(a,g,b,c,form);
fprintf(1,'Done.\n');
figure(design_step); clf
set(design_step,'position',fig_pos{design_step});
DocumentNTF(ABCD,osr,f0);
drawnow;

% Time-domain simulations
fprintf(1,'Doing time-domain simulations... ');
design_step = design_step+1;
figure(design_step); clf;
set(design_step,'position',fig_pos{design_step});
% Example spectrum
subplot('position', [0.05,0.1,0.6 0.8]);
PlotExampleSpectrum(ntf,M,osr,f0);
title('Example Spectrum');
% SQNR plot
subplot('position',[.74 .18 .25 .65]);
if nlev==2
    [snr_pred,amp_pred] = predictSNR(ntf,osr);
    plot(amp_pred,snr_pred,'-');
    hold on;
end
[snr,amp] = simulateSNR(ntf,osr,[],f0,nlev);
fprintf(1,'Done.\n');
plot(amp,snr,'og');
figureMagic([-100 0], 10, 2, [0 100], 10, 2, [7 3], 'Time-Domain Simulations');
xlabel('Input Level (dBFS)');
ylabel('SQNR (dB)');
[peak_snr,peak_amp] = peakSNR(snr,amp);
msg = sprintf('peak SQNR = %4.1fdB  \n@ amp=%4.1fdB  ',peak_snr,peak_amp);
text(peak_amp-10,peak_snr,msg,'hor','right', 'vertical','middle');
msg = sprintf('OSR=%d ',osr);
text(0,5,msg,'hor','right');
title('SQNR Plot');
drawnow;

% Dynamic range scaling
fprintf(1,'Doing dynamic range scaling... ');
design_step = design_step+1;
figure(design_step); clf;
set(design_step,'position',fig_pos{design_step});
ABCD0 = ABCD;
[ABCD umax] = scaleABCD(ABCD0,nlev,f0);
[a,g,b,c] = mapABCD(ABCD,form);
fprintf(1,'Done.\n');
fprintf(1,'Verifying dynamic range scaling... ');
u = linspace(0,0.95*umax,30);
N = 1e4; N0 = 50;
test_tone = cos(2*pi*f0*[0:N-1]);
test_tone(1:N0) = test_tone(1:N0) .* (0.5-0.5*cos(2*pi/N0*[0:N0-1]));
maxima = zeros(order,length(u));
for i = 1:length(u)
    ui = u(i);
    [v,xn,xmax] = simulateDSM( ui*test_tone, ABCD, nlev );
    maxima(:,i) = xmax(:);
    if any(xmax>1e2) 
    	fprintf(1,'Warning, umax from scaleABCD was too high.\n');
    	umax = ui;
    	u = u(1:i);
    	maxima = maxima(:,1:i);
    	break;
    end
end
fprintf(1,'Done.\n');
colors = hsv(order);
cmd = 'legend( handles';
handles = [];
for i = 1:order
    handles(i) = plot(u,maxima(i,:),'--','color',colors(i,:));
    if i==1
    	hold on;
    end
        cmd = [ cmd ',''State ' num2str(i) '''' ];
        plot(u,maxima(i,:),'o','color',colors(i,:));
end
hold off; grid on;
text(umax/2,0.05,'DC input','Hor','cen','Vert','mid');
figureMagic([ 0 umax],[],[], [0 1],0.1,2, [3 3],'State Maxima');
set(gca,'position', [0.1 0.07, 0.85, 0.85]);
cmd = [cmd ',4);'];
eval(cmd);

% The next step would be to size capacitors for each integrator state based 
% on noise (kT/C) limits and area.

% Simulations modelling the effects of finite op-amp gain and capacitor 
% errors should also be performed.

adc.order = order;
adc.osr = osr;
adc.nlev = nlev;
adc.f0 = f0;
adc.ntf = ntf;
adc.ABCD = ABCD;
adc.umax = umax;
adc.peak_snr = peak_snr;
adc.form = form;
adc.coefficients.a = a;
adc.coefficients.g = g;
adc.coefficients.b = b;
adc.coefficients.c = c;
