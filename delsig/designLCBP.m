function [param,H,L0,ABCD,x] = designLCBP(n,OSR,opt,Hinf,f0,t,form,x0,dbg)
%[param,H,L0,ABCD,x]=designLCBP(n=3,OSR=64,opt=2,Hinf=1.6,f0=1/4,t=[0 1],form='FB',x0,dbg=0)
%Design an LC modulator of order 2*n with center frequency f0.
%The feedback waveform is a rectangular pulse from t(1) to t(2).
%Arguments
% n	The number of LC tanks.
% OSR	The oversampling ratio.
% opt	A flag for NTF zero optimization. See synthesizeNTF.
% Hinf	The target out-of-band gain of the NTF.
% f0	The center frequency of the modulator, normalized to fs.
% t	A 2-element vector containing the (normalized) feedback pulse edge times.
% form	The modulator topology. One of 'FB', 'FF', 'FFR' or 'GEN'.
% x0	An initial guess at the parameter vector for the optimizer.
% dbg	A flag for enabling debugging/progress report output.
%
%Output
% param	A struct containing the (n, OSR, Hinf, f0, t and form) arguments
%   plus the following fields:
%   l   The inductances in each tank.
%   c	The capacitances in each tank.
%   gu	The transconductance from the u input to each of the tanks.
%	The final term is the voltage gain from u to the comparator input.
%	(1 by n+1)
%   gv	The transconductance from the v output to each of the tanks. (1 by n)
%   gw	The inter-stage transconductances. (1 by n-1)
%   gx	The gains from the top of each tank resistor to the comparator. (1 by n)
%   rx	The resistances inserted between the output of each interstage 
%	transconductance and top of each tank. Note that r(1) is not used.
%	(1 by n)
%
% H	The noise transfer function of the modulator, taking into account
%	the post-design quantizer gain.
% L0	A description of the open-loop TF from u to v, for use with evalTFP().
% ABCD	The state-space representation of the discrete-time equivalent.

% Handle the input arguments
parameters = {'n';'OSR';'opt';'Hinf';'f0';'t';'form';'x0';'dbg'};
defaults = { 3 64 2 1.6 0.25 [0 1] 'FB' NaN 0 };
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end
param.n = n;	param.OSR = OSR;	param.Hinf = Hinf;
param.f0 = f0;	param.t = t;		param.form = form;

% Compute the resonant frequencies for the LC tanks 
% and the L and C values.
if opt
    w = 2*pi*(f0 + 0.25/OSR*ds_optzeros(n,opt)');
else
    w = 2*pi*f0*ones(1,n);
end
param.l = (1./w);	param.c = (1./w);	% Impedance scaling is possible.

% Find values for the parameters that yield a stable system.
% ?? Could some control theory help here ??
if isnan(x0)
    x = [2:n+1];
elseif length(x0) ~= n
    error('Initial guess (x0) has the wrong number of parameters.');
else
    x = x0;
end

if exist('fmincon','file')==2	% Latest version of Optimization Toolbox
    [param,H,L0,ABCD,x] = designLCBP6(n,OSR,opt,Hinf,f0,t,form,x0,dbg);
elseif exist('constr','file')==2	% Older version of Optimization Toolbox
    if dbg
	fprintf(1, 'Iterations for initial pole placement:\n');
    end
    options = foptions;
    options(2) = 0.001;	% x tolerance
    options(3) = 1;	% f tolerance
    options(14) = 1000;	% Maximum number of iterations
    x = constr('LCObj1',x,options,[],[],[],param,0.97,dbg);
    H = LCoptparam2tf(x,param);
    rmax = max(abs(H.p{:}));
    if rmax>0.97	% Failure to converge!
	fprintf(2,'Unable to find parameter values which stabilize the system.\n');
	return
    end
    % Do the main optimization for the feedback parameters.
    if dbg>1
	options(1) = 1; 	% Extra debugging info
    end
    if dbg
	fprintf(1, '\nParameter Optimization:\n');
    end
    options(2)  = 1e-4;	% x tolerance
    options(3)  = 0.5;	% f tolerance
    options(4)  = 0.01;	% constraint tolerance
    %options(13) = 1;	% Number of equality constraints
    options(14) = 1000;	% Maximum number of iterations.
    options(16) = 1e-4;	% Minimum dx for finite difference gradients
    options(17) = 1e-2;	% Maximum dx for finite difference gradients
    options(18) = .01;	% Step length. (This doesn't prevent large steps!)
    % I have used minimax() and constr() below
    [x,options] = minimax('LCObj',x,options,[],[],[],param,dbg);

    [H L0 ABCD param] = LCoptparam2tf(x,param);

    % Uncomment the following line to take into account the quantizer gain 
    % before returning results and doing the plots
    %[H L0 ABCD k] = LCparam2tf(param,0);
    % Uncomment the following lines to yield parameters corresponding to k=1
    %param.gu = k*param.gu; param.gv = k*param.gv;
    %ABCD(2*n+1,:) = ABCD(2*n+1,:)/k;
    %ABCD(:, 2*n+[1 2]) = ABCD(:, 2*n+[1 2]) *k;

    % Now fix up the inputs for 0dB gain at f0.
    gain_f0 = abs(evalTFP(L0,H,f0)); 
    param.gu = param.gu/gain_f0; L0.k = L0.k/gain_f0;
    ABCD(:,2*n+1) = ABCD(:,2*n+1)/gain_f0;
    if strcmp(form,'FF') | strcmp(form,'FFR')
	param.gu(n+1) =1;		% For the flattest STF
    end

    if dbg
	LCplotTF(H,L0,param);
    end
else
    fprintf(1, 'Sorry, designLCBP needs the optimization toolbox.\n');
end

return
