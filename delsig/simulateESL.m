function [sv,sx,sigma_se,max_sx,max_sy] = simulateESL(v,mtf,M,dw,sx0)
%[sv,sx,sigma_se,max_sx,max_sy] = simulateESL(v,mtf,M=16,dw=[1..],sx0=[0..])
%Simulate the Element Selection Logic for a multi-element DAC.
%Outputs:
% sv	is an MxN matrix whose columns are the selection vectors.
% sx	is an orderxM matrix containing the final state of the ESL.
% sigma_se	is the rms value of the selection error.
% max_sx	is the maximum absolute value of the state for all modulators.
% max_sy	is the maximum absolute value of the input to the VQ.
%Inputs:
% v	is a vector of the digital (positive integer) input values.
% mtf	is the mismatch-shaping transfer function, given in zero-pole form.
% M	is the number of elements.
% sx0	is a matrix whose columns are the initial states of the ESL.
% dw	is a vector of dac element weights
%
% For a brief description of the theory of mismatch-shaping DACs, see
% R. Schreier and B. Zhang "Noise-shaped multibit D/A convertor employing
% unit elements," Electronics Letters, vol. 31, no. 20, pp. 1712-1713,
% Sept. 28 1995.
%

% Handle the input arguments
if nargin < 2
    fprintf(1,'simulateESL needs at least two input arguments\n');
    return
end
parameters = {'v' 'mtf' 'M' 'dw' 'sx0'};
defaults = { NaN NaN 16 NaN NaN };
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
    if strcmp(class(mtf),'zpk')
    end
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
    sx0 = zeros(order,M);
end
if isnan(dw)
    dw = ones(M,1);
end

% B/A = MTF-1
num = poly(mtf.zeros);
den = poly(mtf.poles);
A = real(-den(2:order+1));
B = real(num(2:order+1)+A);

N = length(v);
sv = zeros(M,N);

sx = sx0;
max_sx = max(max(abs(sx)));
max_sy = 0;
sum_se2 = 0;

for i=1:N
    % Compute the sy vector.
    sy = B*sx;
    % Normalize sy for a minimum value of zero.
    sy = sy - min(sy);
    % Pick the elements that have the largest desired_usage (sy) values.
    sv(:,i) = ESLselect(v(i),sy,dw);
    % Compute the selection error.
    se = sv(:,i)' - sy;
    % Compute the new sx matrix
    sxn = A*sx + se;
    sx = [ sxn; sx(1:order-1,:)];
    % Keep track of some statistics.
    sum_se2 = sum_se2 + sum((se-mean(se)).^2);
    max_sx = max([max_sx abs(sxn)]);
    max_sy = max([max_sy abs(sy)]);
end
sigma_se = sqrt(sum_se2/(M*N));
