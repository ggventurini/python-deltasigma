ABCD = [ 1 0 0 0 1 -1 0
         1 1 0 0 0 -2 0
         0 1 1 0 0 0 -1
         0 0 1 1 0 0 -2
         0 1 0 0 0 0 0
         0 0 0 1 0 0 0 ];
nlev = [9 9]; 
[ntf stf] = calculateTF(ABCD, [1 1]);
ncf1 = -ntf(2,1); % 2*(z-0.5)^2/z^4
ncf2 = ntf(1,1); % (z-1)^2/z^2
% stf_eff = stf(1,1)*ncf1 + stf(2,1)*ncf2;
stf_eff= cancelPZ( zpk(tf(stf(1,1)*ncf1) + tf(stf(2,1)*ncf2)) );
ntf_eff = ntf(2,2)*ncf2; % (z-1)^4/z^4

osr = 16;
M = nlev-1;
N = 2^12;
f_bin = round(0.01*N);
t = 0:N-1;
wdw = ds_hann(N)/(M(1)*N/4);

u = undbv(-3)*M(1)*sin( 2*pi*f_bin/N*t);
[v xn xmax ]= simulateDSM(u,ABCD,nlev);
V1 = fft(v(1,:).*wdw);
sqnr1 = calculateSNR(V1(1:N/(2*osr)+1), f_bin);
[freq spec1] = logsmooth(V1,f_bin);

vf = ds_filt([ncf1 ncf2],v);
Vf = fft(vf.*wdw);
sqnr2 = calculateSNR(Vf(1:N/(2*osr)+1), f_bin);
[freq spec2] = logsmooth(Vf,f_bin);

figure(1); clf;
h = semilogx(freq,spec1, 'b');
set(h,'LineWidth',2);
hold on;
h = semilogx(freq,spec2, 'm');
set(h,'LineWidth',2);
msg = sprintf('@ A = %.0f dBFS\n& OSR = %d,\nSQNR1 = %.0f dB\nSQNR = %.0f dB', dbv(Vf(f_bin+1)), osr, sqnr1, sqnr2 );
h = text(0.013,0,msg, 'Vert','Top', 'FontSize',20);
text(0.5,-140,sprintf('NBW = %.1E',1.5/N), 'Hor','right', 'Vert','bot');
figureMagic([1e-3 0.5],[],[],[-140 5],10,2);
plot([1e-3 0.5/osr], -140*[1 1], 'k', 'LineWidth',4);
xlabel('Normalized Frequency');
ylabel('dBFS/NBW');

return;

figure(1); clf
t_plot = 100;
stairs( t(1:t_plot), u(1:t_plot), 'm');
hold on;
stairs( t(1:t_plot), v(1:t_plot), 'b');

figure(2); clf
window = ds_hann(N)/(N/4);
V = fft(v.*window);
semilogx( dbv(V(2:N/2+1)) );
grid on;

