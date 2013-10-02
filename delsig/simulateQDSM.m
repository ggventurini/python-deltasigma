function [v,xn,xmax,y] = simulateQDSM(u,arg2,nlev,x0)
%[v,xn,xmax,y] = simulateQDSM(u,ABCD,nlev=2,x0=0)
% or
%[v,xn,xmax,y] = simulateQDSM(u,ntf,nlev=2,x0=0)
%
%Compute the output of a quadrature delta-sigma modulator with input u,
%a structure described by ABCD, an initial state x0 (default zero) and 
%a quantizer whose number of levels is specified by nlev.
%For multiple quantizers, make nlev a column vector;
%for complex quantization to a diamond lattice, multiply nlev by 1i.
% size(u) = [nu N], size(nlev) = [nq 1], size(ABCD) = [order+nq order+nq+nu]
%
%Alternatively, the modulator may be described by an NTF.
%The NTF is zpk object. (The STF is assumed to be 1.)

if nargin<2
    fprintf(1,'Error. %s needs at least two arguments.\n', mfilename);
    return
end

% Handle the input arguments
parameters = {'u','arg2','nlev','x0'};
defaults = [ NaN NaN 2 NaN ];
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
	eval([parameter '=defaults(i);'])  
    end
end
nu = size(u,1);
nq = length(nlev);
if isobject(arg2) & strcmp(class(arg2),'zpk')
    ntf.k = arg2.k;
    ntf.zeros = arg2.z{:};
    ntf.poles = arg2.p{:};
    form = 2;
    order = length(ntf.zeros);
elseif isstruct(arg2)
    if any(strcmp(fieldnames(arg2),'zeros'))
	ntf.zeros = arg2.zeros;
    else
	error('No zeros field in the NTF.')
    end
    if any(strcmp(fieldnames(arg2),'poles'))
	ntf.poles = arg2.poles;
    else
	error('No poles field in the NTF.')
    end
    form = 2;
    order = length(ntf.zeros);
elseif isnumeric(arg2)
    if size(arg2,2) > 2 & size(arg2,2)==nu+size(arg2,1)	% ABCD dimesions OK
	form = 1;
	ABCD = arg2;
	order = size(ABCD,1)-nq;
    else
	fprintf(1,'The ABCD argument does not have proper dimensions.\n');
	if size(arg2,2) == 2		% Probably old (ver. 2) ntf form
	    fprintf(1,'You appear to be using the old-style form of NTF specification.\n Automatic converstion to the new form will be done for this release only.\n');
	    ntf.zeros = arg2(:,1);
	    ntf.poles = arg2(:,2);
	    form = 2;
	    order = length(ntf.zeros);
	else
	    error('Exiting simulateQDSM.')
	end
    end
else
    error('The second argument is neither an ABCD matrix nor an NTF.\n');
end
if isnan(x0)
    x0 = zeros(order,1);
end

if form==1
    [A B C D] = partitionABCD(ABCD, nq+nu);
    D1 = D(:,1:nu);
else
    % Create a FF realization of 1-1/H. 
    % Note that MATLAB's zp2ss and canon functions don't work for complex TFs.
    A = zeros(order);
    B2 = [1; zeros(order-1,1)];
    diag = 1:order+1:order*order;
    A(diag) = ntf.zeros;
    subdiag = 2:order+1:order*order;
    A(subdiag) = 1;
    % Compute C st C*inv(zI-A)*B = 1-1/H(z);
    w = 2*pi*rand(1,2*order);
    desired = 1-1./evalTF(ntf,exp(1i*w));
    % suppress warnings about complex TFs
    s = warning;
    warning('off');
    sys = ss(A,B2,eye(order),zeros(order,1),1);
    warning(s);
    sysresp = reshape(freqresp(sys,w),order,length(w));
    C = desired/sysresp;
    % !!!! Assume stf=1
    B1 = -B2;
    B = [B1 B2];
    D1 = 1;
end

N = length(u);
v = zeros(nq,N);
y = zeros(nq,N);
if nargout > 1	% Need to store the state information
    xn = zeros(order,N); 
end
if nargout > 2	% Need to keep track of the state maxima
    xmax = abs(x0);
end

for i=1:N
    y(:,i) = C*x0 + D1*u(:,i);
    v(:,i) = ds_qquantize(y(:,i),nlev);
    x0 = A * x0 + B * [u(:,i);v(:,i)];
    if nargout > 1	% Save the next state
	xn(:,i) = x0;
    end
    if nargout > 2	% Keep track of the state maxima
	xmax = max(abs(x0),xmax);
    end
end
return

function v = ds_qquantize(y,n)
%v = ds_qquantize(y,n=2)
%Quadrature quantization
if nargin<2
    n=2;
end
if isreal(n)
    v = ds_quantize(real(y),n) + 1j*ds_quantize(imag(y),n);
else
    ytmp = [real(y)+imag(y) real(y)-imag(y)];
    v = ds_quantize(ytmp,abs(n))*[1+1i;1-1i]/2;
end
return

