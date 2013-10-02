function [f,p] = bplogsmooth(X,tbin,f0)
% Smooth the fft, X, and convert it to dB.
% Use 8 bins from the bin corresponding to f0 to tbin and again as far. 
% Thereafter increase bin sizes by a factor of 1.1, staying less than 2^10.
% For tbin, group the bins together. TBIN IS ASSUMED TO BE IN THE UPPER SIDEBAND!
% Use this for nice double-sided log-log plots
N=length(X);
n=8;

bin0 = round(f0*N);
bin1 = rem(tbin-bin0,n)+bin0;
bind = bin1-bin0;
usb1 = [bin1:n:tbin tbin+3:8:tbin+bind ];
m = usb1(length(usb1))+n;
while m+n/2 < N/2
   usb1 = [usb1 m];
   n = min(n*1.1,2^10);
   m = m+n;
end
usb2 = [usb1(2:length(usb1))-1 N/2];

n=8;
lsb2 = [bin1:-n:bin1-2*bind]-1;
m = lsb2(length(lsb2))-n;
while m-n/2 > 1
   lsb2 = [lsb2 m];
   n = min(n*1.1,2^10);
   m = m-n;
end
lsb1 = [lsb2(2:length(lsb2))+1 1];

startbin=[lsb1(length(lsb1):-1:1) usb1];
stopbin=[lsb2(length(lsb2):-1:1) usb2];

f = ((startbin+stopbin)/2 -1)/N -f0;
p = zeros(size(f));
for i=1:length(f)
    p(i) = dbp(norm(X(startbin(i):stopbin(i)))^2 / ...
              (stopbin(i)-startbin(i)+1));
end
