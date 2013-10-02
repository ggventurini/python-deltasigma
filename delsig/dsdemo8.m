% Demonstrate the designLCBP function

clc;
fprintf(1,'\t\t\t Continuous-Time LC Modulator Design\n\n');
%WARNING: designLCBP is very ill-behaved due to deficiencies
% in MATLAB's constr() and minimax(). These functions frequently
% take a viable set of parameters as a starting point and turn
% them into an unstable modulator. I should provide a more robust
% implementation...
echo on;
n = 3;
OSR = 64;
opt = 2;
Hinf = 1.7;
f0 = 1/16;
t = [0.5 1];
form = 'FB';
dbg = 1;
[param H L0 ABCD] = designLCBP(n,OSR,opt,Hinf,f0,t,form,[],dbg);
echo off;
