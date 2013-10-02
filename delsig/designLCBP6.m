function [param,H,L0,ABCD,x] = designLCBP6(n,OSR,opt,Hinf,f0,t,form,x0,dbg)
% Modified designLCBP for use with latest optimization toolbox

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
    if dbg
	fprintf(1, 'Iterations for initial pole placement:\n');
    end
    options = optimset('TolX',0.001, 'TolFun',1, 'MaxIter',1000 );
    options = optimset(options,'LargeScale','off');
    x = fmincon(@LCObj1a,x,[],[],[],[],[],[],@LCObj1b,options,param,0.97,dbg);
    H = LCoptparam2tf(x,param);
    rmax = max(abs(H.p{:}));
    if rmax>0.97	% Failure to converge!
	fprintf(2,'Unable to find parameter values which stabilize the system.\n');
	return
    end
    % Do the main optimization for the feedback parameters.
    if dbg
	fprintf(1, '\nParameter Optimization:\n');
    end
    options = optimset('TolX',1e-4, 'TolFun',0.5, 'TolCon', 0.01);
    options = optimset(options,'MaxIter',1000 );
    options = optimset(options,'DiffMaxChange',1e-2);
    options = optimset(options,'DiffMinChange',1e-4);
    options = optimset(options,'LargeScale','off');
    if dbg>1
	options = optimset(options,'Display','iter'); 	% Extra debugging info
    end
    x = fmincon(@LCObj2a,x,[],[],[],[],[],[],@LCObj2b,options,param,0.97,dbg);
else
    fprintf(1, 'Sorry, designLCBP needs the optimization toolbox.\n');
    return
end

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
return


function f = LCObj1a(x,param,max_radius,dbg)
% The objective function for the initial optimization
% process used to put the roots of the denominator of the LCBP NTF inside 
% the unit circle.
f = 1;			% No actual objective
return


function [C,Ceq] = LCObj1b(x,param,max_radius,dbg)
% The constraints function for the initial optimization
% process used to put the roots of the denominator of the LCBP NTF inside 
% the unit circle.
H = LCoptparam2tf(x,param);
objective = 1;			% No actual objective
rmax = max(abs(H.p{:}));
C = rmax - max_radius;
Ceq = [];
if dbg
    fprintf(1,'x = [ ');
    fprintf(1, '%.4f ',x);
    fprintf(1,']\n');
    fprintf(1,'rmax = %f\n\n',rmax);
end
return


function objective = LCObj2a(x,param,rmax,dbg)
% The objective function 
if dbg
    fprintf(1, 'x = [');
    fprintf(1, '%8.5f ', x);
    fprintf(1, ']\n');
end
[H L0 ABCD] = LCoptparam2tf(x,param);
% The objective is the in-band noise, in dB
f1 = param.f0 - 0.25/param.OSR;
f2 = param.f0 + 0.25/param.OSR;
objective = dbv(rmsGain(H,f1,f2)/sqrt(3*param.OSR));
fprintf(1, 'N0^2 = %.1fdB\n', objective);
return

function [C, Ceq] = LCObj2b(x,param,rmax,dbg)
% The constraints: r-rmax, Hinf-Hinfdesired, quantizer input,
% and for the FB form, |stf|-1.1.
[H L0 ABCD] = LCoptparam2tf(x,param);
Ceq = [];

% The constraints are on the maximum radius of the poles of the NTF
% the infinity-norm of the NTF, the peaking of the STF (for the 'FB' form),
% and the maximum input to the quantizer.
max_radius = max(abs(H.p{:}));
H_inf = infnorm(H);
stf0 = abs( evalTFP(L0, H, param.f0) );
[tmp1,tmp2,tmp3,y] = simulateDSM(0.5/stf0*sin(2*pi*param.f0*[0:1000]),ABCD);
ymax = max(abs(y));
if strcmp(param.form,'FB')
    f1 = param.f0 - 0.25/param.OSR;
    f2 = param.f0 + 0.25/param.OSR;
    stf1 = abs( evalTFP(L0, H, f1) )/stf0;
    stf2 = abs( evalTFP(L0, H, f2) )/stf0;
    C = [100*(max_radius-rmax) H_inf-param.Hinf stf1-1.01 stf2-1.01 (ymax-5)/10];
    %C = [100*(max_radius-rmax) H_inf-param.Hinf stf1-1.01 stf2-1.01 log(ymax/5)];
else
    C = [100*(max_radius-rmax) H_inf-param.Hinf (ymax-5)/10]; 
end

if dbg
	fprintf(1, 'constraints =');
	fprintf(1, ' %8.5f', C);
	fprintf(1, '\n');
	fprintf(1, 'rmax=%5.3f, Hinf = %4.2f, ymax = %5.3f\n\n', ...
       max_radius, H_inf, ymax );
end
return
