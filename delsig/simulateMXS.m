function [sv,sigma_se,sigma_xe,x] = simulateMXS(v,M,mtf,xtf,alpha,xu,x0,sv0)
%[sv,sigma_se,max_sx,max_sy,xx] = simulateMXS(v,M,mtf,xtf,alpha,xu,xx0,sv0)
%Simulate simultaneous Mismatch- and Transition-Error-Shaping logic for a multi-element DAC.
%Inputs:
% v		digital input sequence (-M:2:M)
% M     number of DAC elements
% mtf	mismatch-shaping transfer function, given in zero-pole form
% xtf	transition-error shaping transfer function, given in zero-pole form
% alpha weighting of the mismatch-shaping loop. 
%       The xe loop is weighted by 1-alpha.
% xu    (1+xu)/2 is the target fraction of elements that should stay the same
% x0	initial state of the ESL (order x M)
% sv0   previous selection vector (Mx1)
%Outputs:
% sv	MxN matrix whose columns are the selection vectors
% x     orderxM matrix containing the final state of the ESL
% sigma_se	rms value of the selection error in the mismatch-shaping loop
% sigma_xe	rms value of the selection error in the transition-error shaping loop.

% Argument checking and default-setting 
ArgumentsAndDefaults = {
  'v'      NaN
  'M'      16
  'mtf'    zpk(1,0,1,1)
  'xtf'    zpk(1,0,1,1)
  'alpha'  0.5
  'xu'     0.5
  'x0'     []
  'sv0'    []
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

order_m = length(mtf.p{1});
order_x = length(xtf.p{1});
order = order_m + order_x;
if isempty(x0)
    x0 = zeros(order,M);
end
if isempty(sv0)
    sv0 = ds_therm(mod(M,2),M);
end
sv0 = sv0(:);  % Ensure sv0 is a column vector

% Bm/Am = MTF-1
num = poly(mtf.z{1});
den = poly(mtf.p{1});
Am = real(-den(2:order_m+1));
Bm = real(num(2:order_m+1)+Am);
% Bx/Ax = XTF-1
num = poly(xtf.z{1});
den = poly(xtf.p{1});
Ax = real(-den(2:order_x+1));
Bx = real(num(2:order_x+1)+Ax);

N = length(v);
sv = zeros(M,N);

sx = x0(1:order_m,:);
xx = x0(order_m+1:end,:);
sum_se2 = 0;
sum_xe2 = 0;

for i = 1:N
    % Compute the xy and sy vectors.
    sy = Bm*sx;
    sy = sy - min(sy);
    xy = Bx*xx + xu;
    yy = alpha*sy + (1-alpha)*xy.*sv0';
    % Pick the elements that have the largest desired_usage (sy) values.
    sv_i = selectElement(v(i),yy);
    sv(:,i) = sv_i;
    % Compute the selection error.
    se = sv_i' - sy;
    xe = (sv_i.*sv0)' - xy;
    % Compute the new state
    sxn = Am*sx + se;
    sx = [ sxn; sx(1:end-1,:)];
    xxn = Ax*xx + xe;
    xx = [ xxn; xx(1:end-1,:)];
    % Keep track of some statistics.
    sum_xe2 = sum_xe2 + sum(xe.^2);
    sum_se2 = sum_se2 + sum((se-mean(se)).^2);
    sv0 = sv(:,i); 
end
sigma_se = sqrt(sum_se2/(M*N));
sigma_xe = sqrt(sum_xe2/(M*N));
x = [sx; xx];
