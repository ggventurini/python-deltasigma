function [objective,constraints] = LCObj(x,param,dbg)
% LCObj	Compute the objective function 
% 	and the constraints: r-rmax, Hinf-Hinfdesired, |stf|-1.
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

% The constraints are on the maximum radius of the poles of the NTF
% the infinity-norm of the NTF, the peaking of the STF (for the 'FB' form),
% and maximum the input to the quantizer.
max_radius = max(abs(H.p{:}));
H_inf = infnorm(H);
stf0 = abs( evalTFP(L0, H, param.f0) );
[tmp1,tmp2,tmp3,y] = simulateDSM(0.5/stf0*sin(2*pi*param.f0*[0:1000]),ABCD);
ymax = max(abs(y));
if strcmp(param.form,'FB')
% 	stf1 = abs( evalTFP(L0, H, f1) )/stf0;
% 	stf2 = abs( evalTFP(L0, H, f2) )/stf0;
% 	constraints = [20*(max_radius-0.97) H_inf-param.Hinf stf1-1.01 stf2-1.01 (ymax-5)/10];
    constraints = [100*(max_radius-0.97) H_inf-param.Hinf log(ymax/5)];
else
    constraints = [20*(max_radius-0.97) H_inf-param.Hinf (ymax-5)/10]; 
end

if dbg
	fprintf(1, 'constraints =');
	fprintf(1, ' %8.5f', constraints);
	fprintf(1, '\n');
	fprintf(1, 'N0=%5.1fdB, rmax=%5.3f, Hinf = %4.2f, ymax = %5.3f\n\n', ...
      objective, max_radius, H_inf, ymax );
end
