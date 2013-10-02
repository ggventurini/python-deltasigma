function mod = mod2()
% mod = mod2()
% Populate a structure describing the second-order modulator
A=[1 0; 1 1];
B=[1 -1; 1 -2];
C=[0 1];
D=[0 0];
ABCD=[A B; C D];
mod.ABCD=[A B; C D];
[mod.H mod.G] = calculateTF(mod.ABCD);
