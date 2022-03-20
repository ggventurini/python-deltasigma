function [sv,sx,max_sy] = simulateTSMS(v,M,mtf,d,sx0)
%[sv,sx,max_sy] = simulateTSMS(v,M=16,mtf,d=0,sx0)
% Simulate tree-structured mismatch-shaping.
%Inputs:
% v      A vector of the input values. v is in (-M:2:M)
% M		 The number of elements.
% mtf	 The mismatch-shaping transfer function, given in zero-pole form.
% d 	 A flag indicating that dither is used.
% sx0	 A matrix whose columns are the initial states of the ESL.
%Outputs:
% sv	 An MxN matrix whose columns are the selection vectors.
% sx	 An orderx(M-1) matrix containing the final state of the ESL.
% max_sy The maximum absolute value of the input to the VQ.%
%
%See 
% I. Galton, "Spectral shaping of circuit errors in digital-to-analog 
% converters," IEEE Trans. Circuits Syst. II, pp. 808-817, Oct. 1997.

% Argument checking and default-setting 
ArgumentsAndDefaults = {
 'v'   NaN
 'M'   16
 'mtf' []
 'd'   0
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

log2M = round(log2(M));
if M ~= 2^log2M
    fprintf(1,'Error. M must be a power of 2.');
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
elseif isstruct(mtf)
    if strcmp(class(mtf),'zpk')
    end
    if ~any(strcmp(fieldnames(mtf),'zeros'))
        error('No zeros field in the MTF.')
    end
    if ~any(strcmp(fieldnames(mtf),'poles'))
        error('No poles field in the MTF.')
    end
elseif isempty(mtf)
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
if isempty(sx0)
    sx0 = zeros(order,M-1);
end

% B/A = MTF-1
num = poly(mtf.zeros);
den = poly(mtf.poles);
A = real(-den(2:order+1));
B = real(num(2:order+1)+A);

N = length(v);
sv = zeros(M,N);

sx = sx0;
max_sy = 0;

vi = zeros(2*M-1,1);
for i=1:N
    maxv = M;
    k1 = 1;
    vi(1) = (v(i)+M)/2; % Translate -M:2:M to 0:M
    for k = 0:log2M-1
        % Compute the intermediate signals (Recursion would be easier!)
        for k2 = 2^(k+1):2:2^(k+2)-1
            y = B*sx(:,k1);
            max_sy = max( max_sy, abs(y) );
            % Implement the even/odd control
            if mod(vi(k1), 2) == 0
                s = 2*round(y/2+d*(rand(1,1)-0.5)); % Round to even integer
            else
                s = 2*round((y-1)/2+d*(rand(1,1)-0.5))+1; % Round to odd integer
            end
            s_limit = min( maxv-vi(k1), vi(k1) );
            s = max( min(s,s_limit), -s_limit);
            e = s-y;
            
            sx(:,k1) = [A*sx(:,k1)+e; sx(1:order-1,k1)];

            vi(k2)= 0.5*(vi(k1) + s);
            vi(k2+1) = 0.5*(vi(k1) - s);
            
            k1 = k1+1;
        end
        maxv = 0.5*maxv;
    end
    sv(:,i) = vi(M:(2*M-1));
end
sv = 2*sv - 1; % Translate [0 1] to [-1 1]
