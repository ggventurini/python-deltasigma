function [sv,sx] = simulateSwap(v,M,sx0)
%[sv,sx] = simulateSwap(v,M=16,sx0=[0..])
%Simulate Adams's Butterfly swapping scheme
%Inputs:
% v         A 1xN vector of the digital input values. v is in (-M:2:M).
% M         The number of elements. M must be a power of 2, i.e. M = 2^m.
% sx0	    An (M/2 x m-1) matrix are the initial states of the ESL.
%Outputs:
% sv	    An MxN matrix whose columns are the selection vectors.
% sx	    An (M/2 x m-1) matrix containing the final state of the ESL.
%
%

% Argument checking and default-setting
ArgumentsAndDefaults = {
    'v'   NaN
    'M'   16
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
m = round(log2(M));
if M ~= 2^m
    error(1,'M=%d must be a power of 2.',M);
end
if isempty(sx0)
    sx0 = zeros(M/2,m);
end

sv = ds_therm(v,M);
sx = sx0;

% Populate P, the permutation matrix for each step
P = [0;1];
pp = 2;
for j=1:m-1
    P = [P; P+pp];
    pp = pp*2;
    Pn = (dec2bin(0:pp-1)-'0')*(2.^(0:j))';
    P = [Pn P];
end
P = P+1;    %Matlab indexes from 1

N = length(v);
for n=1:N
    in = sv(:,n);
    for j=1:m
        for i=1:M/2
            if in(2*i-1) == in(2*i)
                out(2*i+[-1 0]) = in(2*i);
            elseif sx(i,j) == 0
                sx(i,j) = 1;
                out(2*i+[-1 0]) = [1 -1];
            else % sx(i,j) == 1
                sx(i,j) = 0;
                out(2*i+[-1 0]) = [-1 1];
            end
        end
        in = out(P(:,j));
    end
    sv(:,n) = out;
end
