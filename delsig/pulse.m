function y = pulse(S,tp,dt,tfinal,nosum) 
% y = pulse(S,tp=[0 1],dt=1,tfinal=10,nosum=0)	Calculate the sampled pulse response 
% of a ct system. tp may be an array of pulse timings, one for each input.
% Outputs
% y	The pulse response
%
% Inputs
% S	An LTI object specifying the system.
% tp	An nx2 array of pulse timings
% dt	The time increment
% tfinal	The time of the last desired sample
% nosum A flag indicating that the responses are not to be summed 

% Handle the input arguments
parameters = {'S','tp','dt','tfinal','nosum'};
defaults = { NaN, [0 1], 1, 10, 0};
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
	eval([parameter '=defaults{i};'])  
    end
end
% Check the arguments
if S.Ts ~= 0
    fprintf(1, 'Error: S must be a cts-time system.\n');
    return
end

% Compute the time increment
dd = 1;
for i=1:prod(size(tp))
    [x di] = rat(tp(i),1e-3);
    dd = lcm(di,dd);
end
[x ddt] = rat(dt,1e-3);
[x df] = rat(tfinal,1e-3);
delta_t = 1 / lcm( dd, lcm(ddt,df) );
delta_t = max(1e-3, delta_t);	% Put a lower limit on delta_t
y1 = step(S,0:delta_t:tfinal);

nd = round(dt/delta_t);
nf = round(tfinal/delta_t);
ndac = size(tp,1);
ni = size(S.b,2);
if rem(ni,ndac)~=0 
    error('The number of inputs must be divisible by the number of dac timings.');
    % This requirement comes from the complex case, where the number of inputs
    % is 2 times the number of dac timings. I think this could be tidied up.
    return
end
nis = ni/ndac; % Number of inputs grouped together with a common DAC timing
	       % (2 for the complex case)
if ~nosum	% Sum the responses due to each input set
    y = zeros(tfinal/dt+1,size(S.c,1),nis);
else
    y = zeros(tfinal/dt+1,size(S.c,1),ni);
end
for i = 1:ndac
    n1 = round(tp(i,1)/delta_t);
    n2 = round(tp(i,2)/delta_t);
    z1 = [n1 size(y1,2) nis];
    z2 = [n2 size(y1,2) nis];
    yy = [zeros(z1); y1(1:nf-n1+1,:,(i-1)*nis+1:i*nis)] ...
       - [zeros(z2); y1(1:nf-n2+1,:,(i-1)*nis+1:i*nis)];
    yy = yy(1:nd:end,:,:);
    if ~nosum	% Sum the responses due to each input set
	y = y + yy;
    else
	y(:,:,i) = yy; 
    end
end

