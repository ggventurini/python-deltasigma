% the file test_snr_data2.mat may be generated with the following:
clear; clc;
order = 3;
osr = 256;
nlev = 2;
f0 = 0.;
Hinf = 1.25;
form = 'CIFB';

ntf = synthesizeNTF(order, osr, 2, Hinf, f0);
[a1 g1 b1 c1] = realizeNTF(ntf, form);
ABCD = stuffABCD(a1, g1, b1, c1, form);
[snr amp] = simulateSNR(ABCD, osr, [], f0, nlev)
% saving ABCD, snr, amp