function [v,xn,xmax,y] = simulateDSM(u,arg2,nlev,x0)
%[v,xn,xmax,y] = simulateDSM(u,ABCD,nlev=2,x0=0)
% or
%[v,xn,xmax,y] = simulateDSM(u,ntf,nlev=2,x0=0)
%
%Compute the output of a general delta-sigma modulator with input u,
%a structure described by ABCD, an initial state x0 (default zero) and 
%a quantizer with a number of levels specified by nlev.
%Multiple quantizers are implied by making nlev an array,
%and multiple inputs are implied by the number of rows in u.
%
%Alternatively, the modulator may be described by an NTF.
%The NTF is zpk object. (The STF is assumed to be 1.)
%The structure that is simulated is the block-diagional structure used by
%zp2ss.m.

fprintf(1,'Warning: You are running the non-mex version of simulateDSM.\n');
fprintf(1,'Please compile the mex version with "mex simulateDSM.c"\n');

if nargin<2
    fprintf(1,'Error. simulateDSM needs at least two arguments.\n');
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
	    error('Exiting simulateDSM.')
	end
    end
else
    error('The second argument is neither an ABCD matrix nor an NTF.\n');
end
if isnan(x0)
    x0 = zeros(order,1);
end

if form==1
    A = ABCD(1:order, 1:order);
    B = ABCD(1:order, order+1:order+nu+nq);
    C = ABCD(order+1:order+nq, 1:order);
    D1= ABCD(order+1:order+nq, order+1:order+nu);
else
    [A,B2,C,D2] = zp2ss(ntf.poles,ntf.zeros,-1);	% A realization of 1/H
    % Transform the realization so that C = [1 0 0 ...]
    Sinv = orth([C' eye(order)])/norm(C); S = inv(Sinv);
    C = C*Sinv;
    if C(1)<0
	S = -S;
	Sinv = -Sinv;
    end
    A = S*A*Sinv; B2 = S*B2; C = [1 zeros(1,order-1)]; % C=C*Sinv; 
    D2 = 0;
    % !!!! Assume stf=1
    B1 = -B2;
    D1 = 1;
    B = [B1 B2];
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
    v(:,i) = ds_quantize(y(:,i),nlev);
    x0 = A * x0 + B * [u(:,i);v(:,i)];
    if nargout > 1	% Save the next state
	xn(:,i) = x0;
    end
    if nargout > 2	% Keep track of the state maxima
	xmax = max(abs(x0),xmax);
    end
end
return

function v = ds_quantize(y,n)
%v = ds_quantize(y,n)
%Quantize y to 
% an odd integer in [-n+1, n-1], if n is even, or
% an even integer in [-n, n], if n is odd.
%
%This definition gives the same step height for both mid-rise
%and mid-tread quantizers.

if rem(n,2)==0	% mid-rise quantizer
    v = 2*floor(0.5*y)+1;
else 		% mid-tread quantizer
    v = 2*floor(0.5*(y+1));
end

% Limit the output
for qi=1:length(n)	% Loop for multiple quantizers
    L = n(qi)-1;
    i = v(qi,:)>L; 
    if any(i)
	v(qi,i) = L;
    end
    i = v(qi,:)<-L;
    if any(i)
	v(qi,i) = -L;
    end
end
