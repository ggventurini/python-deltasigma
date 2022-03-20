function [sv,sx,sigma_se,max_sx,max_sy] = simulateMS(v,M,mtf,d,dw,sx0)
%[sv,sx,sigma_se,max_sx,max_sy] = simulateMS(v,M=16,mtf=zpk(1,0,1,1),d=0,dw=[1..],sx0=[0..])
%Simulate the vector-based Mismatch-Shaping for a multi-element DAC.
%Inputs:
% v	        A vector of the digital input values. v is in (-M:2:M) if dw=[1..]
%	        otherwise v is in [-sum(dw),sum(dw)].
% M	        The number of elements.
% mtf	    The mismatch-shaping transfer function, given in zero-pole form.
% d 	    Dither uniformly distributed in [-d,d] is added to sy.
% dw	    A vector of dac element weights
% sx0	    A matrix whose columns are the initial states of the ESL.
%Outputs:
% sv	    An MxN matrix whose columns are the selection vectors.
% sx	    An orderxM matrix containing the final state of the ESL.
% sigma_se  The rms value of the selection error.
% max_sx    The maximum absolute value of the state for all modulators.
% max_sy    The maximum absolute value of the input to the VQ.
%
% For a brief description of the theory of mismatch-shaping DACs, see
% R. Schreier and B. Zhang "Noise-shaped multibit D/A convertor employing
% unit elements," Electronics Letters, vol. 31, no. 20, pp. 1712-1713,
% Sept. 28 1995.
%
fprintf(1,'Warning: You are running the non-mex version of simulateMS.\n');
fprintf(1,'Please compile the mex version with "mex simulateMS.c"\n');

% Argument checking and default-setting 
ArgumentsAndDefaults = {
 'v'   NaN
 'M'   16
 'mtf' []
 'd'   0
 'dw'  []
 'sx0' []
  };
for i = 1:size(ArgumentsAndDefaults,1)
    parameter = ArgumentsAndDefaults{i,1};
    if i>nargin || eval(['isempty(' parameter ') ']) || ...
      ( eval(['isnumeric(' parameter ') '])  &&  ...
        eval(['length(' parameter ') <= 1']) && ...
        eval(['isnan(' parameter ')']))
        if isnan(ArgumentsAndDefaults{i,2})
            error('%s: Argument %d (%s) is required.',mfilename, i, parameter )
        else
            eval([parameter '= ArgumentsAndDefaults{i,2};'])
        end
    end
end

if isempty(mtf)
    mtf.zeros = 1;
    mtf.poles = 0;
elseif isobject(mtf)
    if ~strcmp(class(mtf),'zpk')
        error('The mtf is an unsupported type of object.')
    end
    tmp = mtf;	clear mtf;
    mtf.form = 'zp';
    mtf.k = tmp.k;
    mtf.zeros = tmp.z{:};
    mtf.poles = tmp.p{:};
elseif isstruct(mtf)
    if ~any(strcmp(fieldnames(mtf),'zeros'))
        error('No zeros field in the MTF.')
    end
    if ~any(strcmp(fieldnames(mtf),'poles'))
        error('No poles field in the MTF.')
    end
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
if isempty(sx0)
    sx0 = zeros(order,M);
end
if isempty(dw)
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
    sv(:,i) = selectElement(v(i), sy + d*(2*rand(1,M)-1), dw);
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
