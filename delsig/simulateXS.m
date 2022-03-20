function [sv,xx,sigma_xe,max_xx,max_sy] = simulateXS(v,M,xtf,xu,xx0,sv0)
%[sv,sxx,sigma_se,max_sx,max_sy] = simulateXS(v,M,xtf,xu,xx0,sv0)
%Simulate Transition-Error-Shaping for a multi-element DAC.
%Inputs:
% v         A vector of the digital input values (from the range -M:2:M).
% M         The number of DAC elements.
% xtf       The transition-error-shaping transfer function, given in zero-pole form.
% xu        (1+xu)/2 is the target fraction of elements that should stay the same
% xx0       An order x M matrix containing the initial state of the ESL.
% sv0       An Mx1 vector containing the previous selection vector
%Outputs:
% sv        An MxN matrix whose columns are the selection vectors.
% xx        An orderxM matrix containing the final state of the ESL.
% sigma_se  The rms value of the selection error.
% max_sx	The maximum absolute value of the state for all modulators.
% max_sy	The maximum absolute value of the input to the VQ.

% Argument checking and default-setting 
ArgumentsAndDefaults = {
 'v'   NaN
 'M'   16
 'xtf' []
 'xu'  0.5
 'xx0' []
 'sv0' []
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

if isobject(xtf) 
    if ~strcmp(class(xtf),'zpk')
        error('The mtf is an unsupported type of object.')
    end
    tmp = xtf;	clear xtf;
    xtf.form = 'zp';
    xtf.k = tmp.k;
    xtf.zeros = tmp.z{:};
    xtf.poles = tmp.p{:};
elseif isstruct(xtf)
    if strcmp(class(xtf),'zpk')
    end
    if ~any(strcmp(fieldnames(xtf),'zeros'))
        error('No zeros field in the MTF.')
    end
    if ~any(strcmp(fieldnames(xtf),'poles'))
        error('No poles field in the MTF.')
    end
elseif isempty(xtf)
    xtf.form = 'zp';
    xtf.k = 1;
    xtf.zeros = 1;
    xtf.poles = 0;
elseif isnumeric(xtf)
    if size(xtf,2) == 2		% Probably old (ver. 2) ntf form
        fprintf(1,'You appear to be using the old-style form of MTF specification.\n Automatic converstion to the new form will be done for this release only.\n');
        tmp.form = 'zp';
        tmp.k = 1;
        tmp.zeros = xtf(:,1);
        tmp.poles = xtf(:,2);
        xtf = tmp;
    end
end
order = length(xtf.poles);
if isempty(xx0)
    xx0 = zeros(order,M);
end
if isempty(sv0)
    sv0 = ds_therm(mod(M,2),M);
end
sv0 = sv0(:);  % Ensure sv0 is a column vector

% B/A = XTF-1
num = poly(xtf.zeros);
den = poly(xtf.poles);
A = real(-den(2:order+1));
B = real(num(2:order+1)+A);

N = length(v);
sv = zeros(M,N);

xx = xx0;
max_xx = max(max(abs(xx)));
max_sy = 0;
sum_xe2 = 0;

for i=1:N
    % Compute the xy and sy vectors.
    xy = B*xx + xu;
    sy = xy .* sv0';
    % Pick the elements that have the largest desired_usage (sy) values.
    sv(:,i) = selectElement(v(i),sy);
    % Compute the selection error.
    xe = (sv(:,i).*sv0)' - xy;
    % Compute the new sx matrix
    sxn = A*xx + xe;
    xx = [ sxn; xx(1:order-1,:)];
    % Keep track of some statistics.
    sum_xe2 = sum_xe2 + sum(xe.^2);
    max_xx = max([max_xx abs(sxn)]);
    max_sy = max([max_sy abs(sy)]);
    sv0 = sv(:,i); 
end
sigma_xe = sqrt(sum_xe2/(M*N));
