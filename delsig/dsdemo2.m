% Demonstrate simulateDSM and simulateSNR
clc
fprintf(1,'\t\t\tDiscrete-Time Simulation\n');
if exist('LiveDemo','var') == 0
    LiveDemo=0;
end
fig1pos1 = [0 630 300 200];
fig1pos2 = [0 372 777 355];
fig1pos3 = [8 562 300 250];
fig1pos4 = [0 200 777 600];
fig2pos1 = [336 562 450 250];
fig2pos2 = [185 372 662 355];
fig3pos1 = [815 565 254 247];
fig3pos2 = [275 365 520 362];

echo on
OSR = 32;
H = synthesizeNTF(5,OSR,1);
N = 8192;
fB = ceil(N/(2*OSR)); ftest=floor(2/3*fB);
u = 0.5*sin(2*pi*ftest/N*[0:N-1]);	% half-scale sine-wave input
v = simulateDSM(u,H); 
echo off;

figure(1); clf;
t = 0:100;
stairs(t, u(t+1),'r');
hold on;
stairs(t,v(t+1),'g');
axis([0 100 -1.2 1.2]);
xlabel('Sample Number');
ylabel('u, v');
set(gcf,'NumberTitle','off'); 
set(gcf,'Name','Modulator Input & Output');
if LiveDemo
    set(1,'position',fig1pos2);
    changeFig(18,2,12);
else
	fprintf(1,'paused\n');
end
pause
if LiveDemo
    set(1,'position',fig1pos1);
    changeFig;
end

f = linspace(0,0.5,N/2+1);
echo on
spec = fft(v.*ds_hann(N))/(N/4);
echo off;
figure(2); clf;
plot( f, dbv(spec(1:N/2+1)), 'b')
figureMagic([0 0.5], 0.05, 2, [-120 0], 20, 1,[],'Output Spectrum');
xlabel('Normalized Frequency')
ylabel('dBFS')
%title('Output Spectrum');
if LiveDemo
    set(2,'position',fig2pos2);
    changeFig(18,2,12);
else
    fprintf(1,'paused\n');
end
pause

echo on
snr = calculateSNR(spec(3:fB+1),ftest-2);
echo off;
text_handle = text(0.05,-10, sprintf('SNR = %4.1fdB @ OSR = %d',snr,OSR),'vert','middle');
if LiveDemo
    set(text_handle,'fontsize',18);
end
if ~LiveDemo
    fprintf(1,'paused\n');
end
pause

echo on
NBW = 1.5/N;
Sqq = 4 * evalTF(H,exp(2i*pi*f)).^2 / 3;
echo off
hold on;
plot( f, dbp(Sqq*NBW), 'm', 'Linewidth', 2 );
text(0.5, -90, sprintf('NBW = %4.1E x f_s ',NBW),'Hor','right');
legend('Simulation','Expected PSD',4);

if LiveDemo
    pause;
    set(2,'position',fig2pos1);
    changeFig;
    legend('Simulation','Expected PSD',4);
else
    fprintf(1,'paused\n');
end
pause;

echo on
[snr_pred,amp_pred] = predictSNR(H,OSR);
[snr,amp] = simulateSNR(H,OSR);
echo off
figure(3); clf;
plot(amp_pred,snr_pred,'-',amp,snr,'og');
figureMagic([-100 0], 10, 1, [0 100], 10, 1,[],'SQNR');
xlabel('Input Level (dBFS)');
ylabel('SQNR (dB)');
%title('SNR curve- theory and simulation');
[pk_snr pk_amp] = peakSNR(snr,amp);
text(-25,85,sprintf('peak SNR = %4.1fdB\n@ OSR = %d\n',pk_snr,OSR),'Hor','right');
if LiveDemo
    set(3,'position',fig3pos2);
    changeFig(18,2,12);
else
	fprintf(1,'paused\n');
end
pause
if LiveDemo
    set(3,'position',fig3pos1);
    changeFig;
end

clc; 
fprintf(1,'\t\t\t Bandpass Modulator\n');
echo on
f0 = 1/8; OSR=64;
H = synthesizeNTF(8,OSR,1,[],f0);
fB = ceil(N/(2*OSR)); ftest=round(f0*N + 1/3*fB);
u = 0.5*sin(2*pi*ftest/N*[0:N-1]);	% half-scale sine-wave input
v = simulateDSM(u,H); 
echo off;

figure(1); clf;
t = 0:100;
stairs(t, u(t+1),'r');
hold on;
stairs(t,v(t+1),'g');
axis([0 100 -1.2 1.2]);
xlabel('Sample Number');
ylabel('u, v');
%title('Modulator Input & Output');
if LiveDemo
    set(1,'position',fig1pos2);
    changeFig(18,2,12);
else
	fprintf(1,'paused\n');
end
pause
if LiveDemo
    set(1,'position',fig1pos1);
    changeFig;
end

echo on
spec = fft(v.*ds_hann(N))/(N/4);
echo off;
figure(2); clf;
plot(linspace(0,0.5,N/2+1), dbv(spec(1:N/2+1)))
axis([0 0.5 -140 0]);
f1 = round((f0-0.25/OSR)*N);
f2 = round((f0+0.25/OSR)*N);
snr = calculateSNR(spec(f1:f2),ftest-f1+1);
text(0.15,-10, sprintf(' SNR = %4.1f dB @ OSR=%d)',snr,OSR), 'vert','middle');
grid on;
xlabel('Normalized Frequency')
ylabel('dBFS/NBW')
%title('Output Spectrum');
if LiveDemo
    set(2,'position',fig2pos2);
    changeFig(18,2,12);
    pause;
else
    fprintf(1,'paused\n'); pause
end

echo on
NBW = 1.5/N;
Sqq = 4 * evalTF(H,exp(2i*pi*f)).^2 / 3;
echo off
hold on;
plot( f, dbp(Sqq*NBW), 'm', 'Linewidth', 2 );
text(0.5, -90, sprintf('NBW=%4.1E x f_s ',NBW),'Hor','right','Vert','mid');
legend('Simulation','Expected PSD',4);
if ~LiveDemo
    fprintf(1,'paused\n');
end
pause
if LiveDemo
    set(2,'position',fig2pos1);
    changeFig;
    legend('Simulation','Expected PSD',4);
    drawnow;
end

echo on
[snr_pred,amp_pred] = predictSNR(H,OSR,[],f0);
[snr,amp] = simulateSNR(H,OSR,[],f0);
echo off
figure(3); clf;
plot(amp_pred,snr_pred,'-b',amp,snr,'og');
figureMagic([-110 0], 10, 1, [0 110], 10, 1,[],'SQNR');
xlabel('Input Level (dBFS)');
ylabel('SQNR (dB)');
%title('SNR curve- theory and simulation');
[pk_snr pk_amp] = peakSNR(snr,amp);
text(-20,95,sprintf('peak SNR = %4.1fdB\n@ OSR = %d\n',pk_snr,OSR),'Hor','right');
if LiveDemo
    set(3,'position',fig3pos2);
    changeFig(18,2,12);
else
	fprintf(1,'paused\n');
end
pause
if LiveDemo
    set(3,'position',fig3pos1);
    changeFig;
end

figure(3); clf;
figure(2); clf;
figure(1); clf;
colors =['m','b'];
Hinf_list=[2 8];

for i=1:2
	Hinf = Hinf_list(i);
	col = colors(i);
	clc;
	fprintf(1,'\t\t\t 15-step 7th-order Lowpass\n');
	fprintf(1,'\t\t\t Hinf = %.0f\n\n',Hinf);
	echo on
	OSR=8; M = 16;
	H = synthesizeNTF(7,OSR,1,Hinf);
	N = 8192;
	fB = ceil(N/(2*OSR)); ftest=floor(2/7*fB);
	u = 0.5*M*sin(2*pi*ftest/N*[0:N-1]);	% half-scale sine-wave input
	v = simulateDSM(u,H,M+1); 
	echo off;
	
	figure(1); subplot(2,1,i); 
	t = 0:100;
	stairs(t, u(t+1),'g');
	hold on;
	stairs(t,v(t+1),col);
	figureMagic([0 100],10,2, [-M M],2,3,[],'Input & Output');
	xlabel('Sample Number');
	ylabel('u, v');
	%title('Modulator Input & Output');
	if LiveDemo
	    set(1,'position',fig1pos4);
	    changeFig(18,2,12);
	else
		fprintf(1,'paused\n');
	end
	pause
	if LiveDemo
	    set(1,'position',fig1pos3);
	    changeFig;
	end
	
	f = linspace(0,0.5,N/2+1);
	echo on
	spec = fft(v.*ds_hann(N))/(M*N/4);
	echo off;
	figure(2); 
	plot(f, dbv(spec(1:N/2+1)),col)
	snr = calculateSNR(spec(3:fB+1),ftest-2);
	text(0.1,10*(i-3), sprintf(' SNR = %4.1fdB @ OSR=%d ',snr,OSR), ...
	  'vert','middle','Color',col);
	figureMagic([0 0.5], 0.05, 2, [-160 0], 20, 1,[],'Output Spectrum');
	xlabel('Normalized Frequency')
	ylabel('dBFS')
	%title('Output Spectrum');
	if LiveDemo
	    set(2,'position',fig2pos2);
	    changeFig(18,2,12);
	    pause;
	end
	
	echo on
	NBW = 1.5/N;
	Sqq = 4 * evalTF(H,exp(2i*pi*f)).^2 / (3*M^2);
	echo off
	hold on;
	plot( f, dbp(Sqq*NBW), 'c', 'Linewidth', 2 );
	if i==1
		text(0.5, -110, sprintf('NBW=%4.1E x f_s ',NBW),'Hor','right');
	end
 	%legend('Simulation','Expected PSD',4);
	pause
	if LiveDemo
	    set(2,'position',fig2pos1);
	    changeFig;
	    % legend('Simulation','Expected PSD',4); 
	else
	    fprintf(1,'paused\n');
	end
	
	echo on
	[snr,amp] = simulateSNR(H,OSR,[],[],M+1);
	echo off;
	figure(3); 
	plot(amp,snr,['o' col],amp,snr,['--' col]); hold on;
	figureMagic([-120 0], 10, 2, [0 120], 10, 2, [], 'SQNR');
	xlabel('Input Level (dBFS)');
	ylabel('SNR (dB)');
	%title('SNR curve');
	[pk_snr pk_amp] = peakSNR(snr,amp);
	text(-13,pk_snr,sprintf('peak SNR = %4.1fdB\n@ OSR = %d',pk_snr,OSR), ...
	'Hor','right','Color',col);
	if LiveDemo
	    set(3,'position',fig3pos2);
	    changeFig(18,2,12);
	    pause
	    set(3,'position',fig3pos1);
	    changeFig;
	end
end

