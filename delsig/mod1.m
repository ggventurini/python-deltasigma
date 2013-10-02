function mod = mod1()
% mod = mod1()
% Populate a structure describing the first-order modulator
A=1;
B=[1 -1];
C=1;
D=[0 0];
mod.ABCD=[A B; C D];
[mod.H mod.G] = calculateTF(mod.ABCD);
