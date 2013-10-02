function y = undbv(x)
% y = undbv(x)	Convert x from dB to a voltage
y = 10.^(x/20);
