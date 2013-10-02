function [H,L0,ABCD,k]=LCparam2tf(param,k)
%[H,L0,ABCD,k]=LCparam2tf(param,k=1)
%Compute the transfer functions (and the Q. gain, k) of an LC modulator.
%
%Arguments
% param	A struct containing the following fields (all normalized):
%   l   The inductances in each tank.
%   c	The capacitances in each tank.
%   gu	The transconductance from the u input to each of the tanks.
%	The final term is the voltage gain to the comparator input.
%	Default = [1 0 .. ]. (1 by n+1)
%   gv	The transconductance from the v output to each of the tanks.
%	(1 by n+1).
%   gw	The inter-stage transconductances. 
%	Default = [1 1 .. ]. (1 by n-1)
%   gx	The gains from the top of each tank resistor to the comparator.
%	Default = [0 0 .. 1]. (1 by n)
%   rx	The resistances inserted between the output of each interstage 
%	transconductance and the top of each tank.
%	Default = [0 ..]. (1 by n)  Note that rx(1) is not used.
%   t	A two-element vector containing the start and end times of the feedback
%	waveform.
%	Default = [0 1].
%   rl	The series resistance of each inductor. 
%	Default = [0 ..]. (1 by n)
%   gc	The parallel conductances of each tank. 
%	Default = [0 ..]. (1 by n)
%
% k	The value to assume for the quantizer gain. If k is omitted or if
%	k is defaulted, k is taken to be 1.
%	If k is 0, a value is computed by the formula 
%	k = mean(abs(y))/mean(y.^2), where y is the quantizer input
%	sequence when the modulator is fed with a small input at f0.
%
%Output
% H	The closed-loop noise transfer function (in z).
% L0	The open-loop tf (in s) from u to v.  G  may be evaluated
%	using evalContSTF(L0,H,f).
% k	The value of the quantizer gain used to compute the above tfs.
% ABCD  The state-space description of the discrete-time system
%
%For the conversion from continuous to discrete time, corrections to the
%standard MATLAB formulae need to be applied in order for the sampled-data
%value to be taken at the END of the time interval.
%This effectively makes t1>0 and requires corrections if t2>=1 and gv(n+1)~=0.
%More complex formulae are needed for t2>=1, so for the moment this code
%requires t2<1.

l = param.l; 
c = param.c;
f0 = mean(1 ./ sqrt(l.*c)) / (2*pi);
n = length(l);
%Handle the input arguments and their defaults
% !! isempty() must be first because when k=[] isnan returns [] and this 
% !! makes the whole expression []
if nargin<2 | isempty(k) | isnan(k) 
    k = 1;
end
if any(strcmp(fieldnames(param),'gu')) & ~isempty(param.gu) & ~isnan(param.gu)	
    gu = param.gu;
else
    gu = [1 zeros(1,n)];
end
gu = padr(gu,n+1,0);	% Pad out with zeros, if needed.
gv = padr(param.gv,n+1,0);	% Pad out with zeros, if needed.
if any(strcmp(fieldnames(param),'gw')) & ~isempty(param.gw) & ~isnan(param.gw)	
    gw = param.gw;
else
    gw = ones(1,n-1);
end
if any(strcmp(fieldnames(param),'gx')) & ~isempty(param.gx) & ~isnan(param.gx)	
    gx = param.gx;
else
    gx = [zeros(1,n-1) 1];
end
if any(strcmp(fieldnames(param),'rx')) & ~isempty(param.rx) & ~isnan(param.rx)	
    rx = param.rx;
    if length(rx)==1;	
	rx = rx(ones(1,n));
    end
else
    rx = [zeros(1,n)];
end
if any(strcmp(fieldnames(param),'t')) & ~isempty(param.t) & ~isnan(param.t)	
    t = param.t;
else
    t = [0 1];
end
if any(strcmp(fieldnames(param),'rl')) & ~isempty(param.rl) & ~isnan(param.rl)	
    rl = param.rl;
    if length(rl)==1;	
	rl = rl(ones(1,n));
    end
else
    rl = [zeros(1,n)];
end
if any(strcmp(fieldnames(param),'gc')) & ~isempty(param.gc) & ~isnan(param.gc)	
    gc = param.gc;
    if length(gc)==1;	
	gc = gc(ones(1,n));
    end
else
    gc = [zeros(1,n)];
end

% Form the state-space description of the modulator
% The states are ordered v1, i1, v2, i2 ...
n2 = 2*n;
ABc = zeros(n2,n2+2); CDc = zeros(1,n2+2);
for i = 1:n
    i2 = 2*i-1;
    % y represents the voltage at the top of the resistor stacked on the tank.
    if i == 1
	y = [ 1 zeros(1,n2+1) ];
    else
	ABc(i2,:) = gw(i-1)*y;
	y = rx(i)*gw(i-1)*y; y(i2) = 1;
    end
    ABc(i2,n2+1:n2+2) = [gu(i) -gv(i)]; 
    ABc(i2:i2+1,i2:i2+1) = [-gc(i) -1; 1 -rl(i)];
    ABc(i2,:) = ABc(i2,:)/c(i); 
    ABc(i2+1,:) = ABc(i2+1,:)/l(i); 
    CDc = CDc + gx(i)*y;
end
CDc = CDc + [zeros(1,n2) gu(n+1) -gv(n+1)];

% THIS IS THE NEW CODE !!
[Ac,Bc,Cc,Dc] = partitionABCD([ABc;CDc], 2);
sys_c = ss( Ac, Bc(:,2), Cc, Dc(2) );
sys = mapCtoD(sys_c,t,f0);
% Augment the input matrix. 
A=sys.a;
B=[padb(Bc(:,1),size(sys.b,1)) sys.b];
C=sys.c;
D=[Dc(1) sys.d];

% Compute L0; use the LC parameters to compute the poles and thereby
% ensure exact agreement with the zeros in H.
%!! THIS IS GOING TO NEED FIXING
s_poles = zeros(1,2*n);
for i=1:n
    s_poles(2*i-1:2*i) = roots([ l(i)*c(i) ( rl(i)*c(i) + gc(i)*l(i) ) 1+rl(i)*gc(i) ]);
end
LF0 = ss(Ac,Bc(:,1),Cc,Dc(1));
L0 = zpk( LF0 );
L0.p =  s_poles;

% Compute H. Use the poles in L0 to compute the zeros in H.
ABCD =[A B; C D];
if k==0		% Estimate k by simulatiing the modulator with a small input
    w0 = mean(1./sqrt(l.*c));
    H = calculateTF(ABCD);
    H.z = exp(s_poles);
    stf0 = abs(evalTFP(L0,H,w0/(2*pi)));
    u = 0.1/stf0*sin(w0*[0:10000]);
    [tmp1,tmp2,tmp3,y] = simulateDSM(u,[A B; C D]);
    k = mean(abs(y))/mean(y.^2);
end
H = calculateTF(ABCD,k);
H.z = exp(s_poles);

% Correct L0k to include the quantizer gain
L0.k = L0.k*k;
