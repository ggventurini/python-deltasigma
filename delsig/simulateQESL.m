function [sv,sx,sigma_se,max_sx,max_sy] = simulateQESL(v,mtf,M,sx0)
%[sv,sx,sigma_se,max_sx,max_sy] = simulateQESL(v,mtf,M=16,sx0=[0..])
%Simulate the Element Selection Logic for a quadrature differential DAC.
%Outputs:
% sv	is an 2*MxN matrix whose columns are the selection vectors.
% sx	is an orderxM matrix containing the final state of the ESL.
% sigma_se	is the rms value of the selection error.
% max_sx	is the maximum absolute value of the state for all modulators.
% max_sy	is the maximum absolute value of the input to the VQ.
%Inputs:
% v	is a vector of the digital input values.
% mtf	is the mismatch-shaping transfer function, given in zero-pole form.
% M	is the number of elements for each channel. 
%       There are 2*M total elements.
% sx0	is a matrix whose columns are the initial states of the ESL.
%

% Handle the input arguments
if nargin < 2
    fprintf(1,'simulateESL needs at least two input arguments\n');
    return
end
parameters = {'v' 'mtf' 'M' 'sx0'};
defaults = { NaN NaN 16 NaN };
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end

if isobject(mtf) 
    if ~strcmp(class(mtf),'zpk')
	error('The mtf is an unsupported type of object.')
    end
    tmp = mtf;	clear mtf;
    mtf.form = 'zp';
    mtf.k = tmp.k;
    mtf.zeros = tmp.z{:};
    mtf.poles = tmp.p{:};
    form = 2;
elseif isstruct(mtf)
    if ~any(strcmp(fieldnames(mtf),'zeros'))
	error('No zeros field in the MTF.')
    end
    if ~any(strcmp(fieldnames(mtf),'poles'))
	error('No poles field in the MTF.')
    end
elseif isnan(mtf)
    mtf.form = 'zp';
    mtf.k = 1;
    mtf.zeros = 1;
    mtf.poles = 0;
elseif isnumeric(mtf)
    if size(mtf,2) == 2		% Probably old (ver. 2) ntf form
	fprintf(1,'You appear to be using the old-style form of MTF specification.\n Automatic converstion to the new form will be done for this release only.\n');
	tmp.form = 'zp';
	tmp.k = 1;
	tmp.zeros = mtf(:,1);
	tmp.poles = mtf(:,2);
	mtf = tmp;
    end
end
order = length(mtf.poles);
if isnan(sx0)
    sx0 = zeros(order,2*M);
end

% B/A = MTF-1
num = poly(mtf.zeros);
den = poly(mtf.poles);
A = -den(2:order+1);
B = num(2:order+1)+A;

N = length(v);
sv = zeros(2*M,N);

sx = sx0;
max_sx = max(max(abs(sx)));
max_sy = 0;
sum_se2 = 0;

for i=1:N
    %fprintf(1,'\r%4d',i);
    % Compute the sy vector.
    sy = B*sx;
    % Bias sy according to the input data, v
    sy = sy - mean(sy) + v(i)/M;
    % Pick the elements that have the largest desired_usage (sy) values.
    sv(:,i) = selectQESL(v(i),sy);
    % Compute the selection error.
    se = sv(:,i).' - sy;
    % Compute the new sx matrix
    sxn = A*sx + se;
    sx = [ sxn; sx(1:order-1,:)];
    % Keep track of some statistics.
    sum_se2 = sum_se2 + norm(se-mean(se))^2;
    max_sx = max([max_sx abs(sxn)]);
    max_sy = max([max_sy abs(sy)]);
end
sigma_se = sqrt(sum_se2/(2*M*N));


