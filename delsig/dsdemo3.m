% Realization and dynamic range scaling
clc
if exist('LiveDemo','var') == 0
    LiveDemo=0;
end
 
fprintf(1,'\t\t\t Modulator realization and scaling\n\n');
echo on;
order=5; R=42; opt=1;
H = synthesizeNTF(order,R,opt);
[a,g,b,c] = realizeNTF(H);
b = [b(1) zeros(1,length(b)-1)];	% Use a single feed-in for the input
echo off

figure(1); clf;
plotPZ(H);
set(1,'name','NTF');
%title('Poles and zeros of a 5th-order modulator')
if LiveDemo
    set(1,'position',[10 307 480 420]);
    changeFig(18,2,12);
	pause
    set(1,'position',[9 526 200 200]);
    changeFig;
else
	drawnow;	
end

fprintf(1,'\nUnscaled modulator\n');
fprintf(1,'   DAC feedback coefficients = ');
for i=1:order
    fprintf(1,' %.6f',a(i));
end
fprintf(1,'\n   resonator feedback coefficients = ');
for i=1:order/2
    fprintf(1,' %.6f',g(i));
end
fprintf(1,'\n');

fprintf(1,'\nCalculate the state maxima.\n');
ABCD = stuffABCD(a,g,b,c);
echo on;
u = linspace(0,0.6,30);
N = 1e4; 
echo off;
T = ones(1,N);
maxima = zeros(order,length(u));
for i = 1:length(u)
    ui = u(i);
    [v,xn,xmax] = simulateDSM( ui(T), ABCD );
    maxima(:,i) = xmax(:);
    if any(xmax>1e2) 
	umax = ui;
	u = u(1:i);
	maxima = maxima(:,1:i);
    	break;
    end
end
figure(2); clf
for i = 1:order
    semilogy(u,maxima(i,:),'o');
    if i==1
	hold on;
    end
    semilogy(u,maxima(i,:),'--');
end
grid on;
xlabel('DC input')
set(gcf,'NumberTitle','off'); 
set(gcf,'Name','Simulated State Maxima');
axis([ 0 0.6 1e-4 10]);
if LiveDemo
    set(2,'position',[238 355 511 372]);
    changeFig(18,2,8);
	pause
else
	fprintf(1,'paused\n');
	pause;
end

clc
fprintf(1,'\nCalculate the scaled coefficients.\n');
echo on;
[ABCDs,umax] = scaleABCD(ABCD,[],[],[],[],[],1e4);
[as,gs,bs,cs] = mapABCD(ABCDs);
echo off;
fprintf(1,'\nScaled modulator, umax=%.2f\n', umax);
fprintf(1,'   DAC feedback coefficients = ');
for i=1:order
    fprintf(1,' %.6f',as(i));
end
fprintf(1,'\n   resonator feedback coefficients = ');
for i=1:order/2
    fprintf(1,' %.6f',gs(i));
end
fprintf(1,'\n   interstage coefficients = ');
for i=1:order
    fprintf(1,' %.6f',cs(i));
end
fprintf(1,'\n   feed-in coefficients = ');
for i=1:order
    fprintf(1,' %.6f',bs(i));
end
fprintf(1,'\n');

fprintf(1,'\nCalculate the state maxima.\n');
echo on;
u = linspace(0,umax,30);
N = 1e4; 
echo off;
T = ones(1,N);
maxima = zeros(order,length(u));
for i = 1:length(u)
    ui = u(i);
    [v,xn,xmax] = simulateDSM( ui(T), ABCDs );
    maxima(:,i) = xmax(:);
    if any(xmax>1e2) 
	umax = ui;
	u = u(1:i);
	maxima = maxima(:,1:i);
    	break;
    end
end

figure(2); clf
for i = 1:order
    semilogy(u,maxima(i,:),'o');
    if i==1
    	hold on;
    end
    semilogy(u,maxima(i,:),'--');
end
grid on;
xlabel('DC input')
axis([ 0 0.6 4e-2 4]);
if LiveDemo
    set(2,'position',[238 355 511 372]);
    changeFig(18,2,8);
	pause
    set(2,'position',[238 527 494 200]);
    changeFig;
end

