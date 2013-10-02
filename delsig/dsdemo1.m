% Demonstrate synthesizeNTF
if exist('LiveDemo','var') == 0
    LiveDemo=0;
end
format compact;
fig1pos1 = [9 630 200 200];
fig1pos2 = [10 407 480 420];
fig2pos1 = [239 630 450 200];
fig2pos2 = [241 341 523 485];

clc
fprintf(1,'\t\tNTF Synthesis-- 5th-order modulator\n\n');
echo on
order = 5;
OSR = 32;
opt = 0;
H = synthesizeNTF(order,OSR,opt);
echo off

figure(1); clf
plotPZ(H);
set(gcf,'MenuBar','none'); 
set(gcf,'NumberTitle','off'); 
set(gcf,'Name','NTF Poles and Zeros');
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

figure(2); clf
f = [linspace(0,0.75/OSR,100) linspace(0.75/OSR,0.5,100)];
z = exp(2i*pi*f);
magH = dbv(evalTF(H,z));
subplot(211);
plot(f,magH);
figureMagic([0 0.5],0.05,2, [-100 10],10,2 );
xlabel('Normalized frequency (1\rightarrow f_s)');
ylabel('dB');
set(gcf,'MenuBar','none'); 
set(gcf,'NumberTitle','off'); 
set(gcf,'Name','NTF Magnitude Response');
fstart = 0.01;
f = linspace(fstart,1.2,200)/(2*OSR); z = exp(2i*pi*f);
magH = dbv(evalTF(H,z));
subplot(212);
semilogx(f*2*OSR,magH);
axis([fstart 1.2 -100 -30]);
grid on
sigma_H = dbv(rmsGain(H,0,0.5/OSR));
hold on;
semilogx([fstart 1], sigma_H*[1 1]);
plot([fstart 1], sigma_H*[1 1],'o');
text( 0.15, sigma_H+5, sprintf('rms gain = %5.0fdB',sigma_H));
xlabel('Normalized frequency (1\rightarrow f_B)');
ylabel('dB')
if LiveDemo
    set(2,'position',fig2pos2);
    changeFig(18,2,12);
else
	fprintf(1,'paused\n');
end
pause
if LiveDemo
    set(2,'position',fig2pos1);
    changeFig;
end

fprintf(1,'\t\t\tOptimized zeros\n\n');
echo on
opt = 1;
H = synthesizeNTF(order,OSR,opt);
echo off

figure(1); clf
plotPZ(H)
set(gcf,'Name','NTF Poles and Optimized Zeros');

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

figure(2)
fmt = '--';
f = [linspace(0,0.75/OSR,100) linspace(0.75/OSR,0.5,100)];
z = exp(2i*pi*f);
magH = dbv(evalTF(H,z));
subplot(211);
hold on;
plot(f,magH,fmt);

f = linspace(fstart,1.2,200)/(2*OSR); z = exp(2i*pi*f);
magH = dbv(evalTF(H,z));
subplot(212);
semilogx(2*f*OSR,magH,fmt);
sigma_H = dbv(rmsGain(H,0,0.5/OSR));
plot([fstart 1], sigma_H*[1 1],fmt);
plot([fstart 1], sigma_H*[1 1],'o');
text( 0.15, sigma_H+5, sprintf('rms gain = %5.0fdB',sigma_H));
if LiveDemo
    set(2,'position',fig2pos2);
    changeFig(18,2,12);
else
	fprintf(1,'paused\n');
end
pause
if LiveDemo
    set(2,'position',fig2pos1);
    changeFig;
end


clc
fprintf(1,'\t\tNTF Synthesis-- Bandpass Modulator\n\n');
echo on
order = 8;
OSR = 64;
opt = 2;
f0 = 0.125;	% fs/8
H = synthesizeNTF(order,OSR,opt,[],f0);
echo off

figure(1); clf
plotPZ(H)
set(gcf,'Name','Bandpass NTF Poles and Zeros');
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

figure(2); clf
f = [ linspace(0,f0-1/(2*OSR),50) linspace(f0-1/(2*OSR),f0+1/(2*OSR),100) linspace(f0+1/(2*OSR),0.5,50)];
z = exp(2i*pi*f);
magH = dbv(evalTF(H,z));
subplot(211);
plot(f,magH);
hold on
G = zpk(zeros(1,order/2),H.p,1,1);
G.k = 1/abs(evalTF(G,exp(2i*pi*f0)));
magG = dbv(evalTF(G,z));
plot(f,magG,'r');
axis([0 0.5 -100 10]);
grid on;
xlabel('Normalized frequency (1 \rightarrow fs)');
ylabel('dB')
set(gcf,'Name','Bandpass NTF/STF Magnitude Response');

f = linspace(f0-0.3/OSR,f0+0.3/OSR); z = exp(2i*pi*f);
magH = dbv(evalTF(H,z));
subplot(212);
plot(2*OSR*(f-f0),magH);
axis([-0.6 0.6 -100 -60]);
grid on
sigma_H = dbv( rmsGain(H,f0-0.25/OSR,f0+0.25/OSR) );
hold on;
semilogx([-0.5 0.5], sigma_H*[1 1]);
plot([-0.5 0.5], sigma_H*[1 1],'o');
text( -0.2, sigma_H+5, sprintf('rms gain = %5.0fdB',sigma_H));
xlabel('Normalized frequency offset');
ylabel('dB')
if LiveDemo
    set(2,'position',fig2pos2);
    changeFig(18,2,12);
end
pause
if LiveDemo
    set(2,'position',fig2pos1);
    changeFig;
end
