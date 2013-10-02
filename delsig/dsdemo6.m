% Demonstrate Saramaki half-band filter design
if exist('LiveDemo','var') == 0
    LiveDemo=0;
end
J = sqrt(-1);
format compact;

clc
fprintf(1,'\t\tHalf-band filter design\n\n');
fprintf(1,'\t\tALPHA VERSION\n\n');

if 0
    % Use one of the canned examples
    [f1 f2] = exampleHBF(2);
else
    fp = 0.9*0.25;
    delta = undbv( -100 );
    figure(1); clf
    set(gcf,'Name','designHBF Iterations');
    [f1,f2,info]=designHBF(fp,delta,1);
end
n1 = length(f1);
n2 = length(f2);
complexity =  size([f1.csd],2) + (2*n1-1)*(n2+size([f2.csd],2)-1);


% interleave the even and odd decimated impulse responses
Nimp = 2^11;
imp = simulateHBF([1 zeros(1,Nimp-1)],f1,f2);

mag = abs(fft(imp));
mag = mag(1:end/2+1);
figure(2); clf;
set(gcf,'Name','designHBF Result');
plot(linspace(0,0.5,length(mag)),dbv(mag));
figureMagic([0 0.5],0.05,2, [-150 3],10,5);

