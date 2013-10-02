function adc = dsexample4
% Example quadrature bandpass modulator 

% Design parameters
order = 4;
osr = 32;
M = 8;
NG = -50;
ING = -10;
f0 = 1/16;
quadrature = 1;
form = 'PFB';

% Derived paramters
nlev = M+1;
z0 = exp(1j*2*pi*f0);
bw = 1/osr;	% Two-sided BW
adc.order = order;
adc.osr = osr;
adc.M = M;
adc.f0 = f0;
adc.quadrature = quadrature;
adc.form = form;

delta = 2;
FullScale = M;

clc
fprintf(1,'\t\t\t%dth-Order Quadrature Example\n\n',order);
fig_pos = { [  10 595  600 215],
            [ 640 595  600 215],
            [  10 345  450 205],
            };

% NTF synthesis and realization
fprintf(1,'Choosing the NTF and realization... ');
design_step = 1;
ntf0 = synthesizeQNTF(order,osr,f0,NG,ING);
% ntf0.z{1} = [z0 z0 z0 ntf0.z{1}(4)];
ABCD = realizeQNTF(ntf0,form,1);
adc.ABCD = ABCD;
[ntf stf] = calculateTF(ABCD);     % Should check that ntf==ntf0
adc.ntf = ntf;
adc.stf = stf;
fprintf(1,'Done.\n');
figure(design_step); clf;
set(design_step,'position',fig_pos{design_step});
DocumentNTF(adc);
drawnow;
%printmif(fullfile('MIF','ntf'), [6 2], 'Helvetica10')

% Time-domain simulations
fprintf(1,'Doing time-domain simulations... ');
design_step = design_step+1;
figure(design_step); clf;
set(design_step,'position',fig_pos{design_step});
% Example spectrum
subplot('position', [0.05,0.1,0.6 0.8]);
PlotExampleSpectrum(ntf,M,osr,f0,quadrature);
title('Example Spectrum');
drawnow;
% SQNR plot
subplot('position',[.74 .18 .25 .65]);
[snr,amp] = simulateSNR(ABCD,osr,[],f0,nlev);
fprintf(1,'Done.\n');
plot(amp,snr,'ob');
hold on;
figureMagic([-100 0], 10, 2, [0 100], 10, 2, [7 3], 'Discrete-Time Simulation');
xlabel('Input Level (dBFS)');
ylabel('SQNR (dB)');

[peak_snr,peak_amp] = peakSNR(snr,amp);
adc.peak_snr = peak_snr;
msg = sprintf('peak SQNR = %4.1fdB  \n@ amp=%4.1fdB  ',peak_snr,peak_amp);
text(peak_amp-10,peak_snr,msg,'hor','right', 'vertical','middle');
msg = sprintf('OSR=%d ',osr);
text(0,5,msg,'hor','right');
title('SQNR Plot');
drawnow;

% Example I/Q mismatch
fprintf(1,'Calculating effect of example I/Q mismatch... ');
design_step = design_step+1;
figure(design_step); clf;
set(design_step,'position',fig_pos{design_step});
ABCDr = mapQtoR(ABCD);
ABCDr(2,end) = 1.01*ABCDr(2,end);   % 0.1% mismatch in first feedback
[H G HI GI] = calculateQTF(ABCDr);
f = ds_freq(osr,f0,quadrature);
w = 2*pi*f;
NTF  = squeeze( freqresp(H,w) );
INTF = squeeze( freqresp(HI,w) );
STF  = squeeze( freqresp(G,w) );
ISTF = squeeze( freqresp(GI,w) );
fprintf(1,'Done.\n');
plot(f,dbv(NTF),'b', 'Linewidth', 2);
hold on;
plot(f,dbv(INTF),'c', 'Linewidth', 2);
plot(f,dbv(STF),'m', 'Linewidth', 2);
plot(f,dbv(ISTF),'r', 'Linewidth', 2);
inband = find( f>=f0-0.5/osr & f<=f0+0.5/osr );
ng = dbp( mean(abs(NTF(inband)).^2) );
plot(f0+0.5*[-1 1]/osr, ng*[1 1], 'k', 'Linewidth',3 );
msg = sprintf('  NG = %.0fdB ', ng);
text(f0+0.5/osr,ng,msg,'Hor','left','Vert','mid');
ing = dbp( mean(abs(INTF(inband)).^2) );
plot(f0+0.5*[-1 1]/osr, ing*[1 1], 'k', 'Linewidth',3 );
msg = sprintf('  ING = %.0fdB ', ing);
text(f0+0.5/osr,ing,msg,'Hor','left','Vert','mid');
irr = min( dbv(STF(inband)) - dbv(ISTF(inband)) );
plot(f0+0.5*[-1 1]/osr, -irr*[1 1], 'k', 'Linewidth',3 );
msg = sprintf('  IRR = %.0fdB ', irr);
text(f0+0.5/osr,-irr,msg,'Hor','left','Vert','mid');
figureMagic([-0.5,0.5],1/16,2, [-80 15],10,2,[],'Example I/Q Mismatch');
xlabel('frequency');
legend('NTF','INTF','STF','ISTF',3);
drawnow;
%printmif(fullfile('MIF','itf'), [6 2], 'Helvetica10')

return
