function ABCD = realizeQNTF(ntf,form,rot,bn)
% ABCD = realizeQNTF(ntf,form='FB',rot=0,bn)
% Convert a quadrature noise transfer function (ntf)  
% into an ABCD matrix for the desired structure.
% Supported structures are
%	FB	Feedback
%	PFB	Parallel feedback
% 	FF	Feedforward  (bn is the coefficient of the aux DAC)
%   PFF Parallel feedforward
%   FBD	FB with delaying quantizer  NOT SUPPORTED YET
% 	FFD	FF with delaying quantizer  NOT SUPPORTED YET
% rot =1 means rotate states to make as many coefficients as possible real 
% 
% The order of the zeros is used when mapping the NTF onto the chosen topology.

% The basic idea is to equate the value of the loop filter at a set of
% points in the z-plane to the values of L1 = 1-1/H at those points.

% Handle the input arguments
parameters = {'ntf';'form';'rot';'bn'};
defaults = {NaN, 'FB', 0,0};
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end

% Code common to all forms
ntf_p = ntf.p{1};
ntf_z = ntf.z{1};
order = length(ntf_p);

A = diag(ntf_z);
A(2:order+1:end) = 1;
if strcmp(form,'PFB') | strcmp(form,'PFF')
    % partition the NTF zeros into in-band and image-band
    fz = angle(ntf_z)/(2*pi);
    f0 = fz(1);
    inband_zeros = abs(fz-f0)<abs(fz+f0);
    n_in = sum(inband_zeros);
    imageband_zeros = ~inband_zeros;
    n_im = sum(imageband_zeros);
    if any(imageband_zeros(n_in+1:end)~=1) 
	error('Please put the image-band zeros at the end of ntf.z');
    end
    if n_im >0
	A(n_in+1,n_in)=0;
    end
end
D = zeros(1,2);

% Find a set of points in the z-plane that are not close to zeros of H
zSet = zeros(order*2,1);
for i = 1:order*2
    z = 2*rand(1,1)-1 + 1i* (2*rand(1,1)-1);
    while any(abs(z-[ntf_z; zSet])<0.1)
        z = 2*rand(1,1)-1 + 1i* (2*rand(1,1)-1);
    end
    zSet(i) = z;
end
% Evaluate L1 = 1-1/H at these points
L1 = 1 - 1./evalTF(ntf,zSet);

switch form
    case 'FB'
	B = zeros(order,2);
	C = [zeros(1,order-1) 1];
	% Compute F = C * inv(zI-A) at each z in zSet
	F = zeros(length(zSet),order);
	I = eye(order);
	for i = 1:length(zSet)
	    F(i,:) = C * inv(zSet(i)*I-A);
	end
	B(:,2) = F\L1;
	if rot==1
	    ABCD = [A B(:,2); C 0];
	    for i = 1:order
		phi = angle(ABCD(i,end));
		ABCD(i,:) = ABCD(i,:)*exp(-1i*phi);
		ABCD(:,i) = ABCD(:,i)*exp(1i*phi);
	    end
	    [A B2 C D2] = partitionABCD(ABCD);
	    B(:,2) = B2;
	end
	B(1,1) = abs(B(1,2));
    case 'PFB'
	B = zeros(order,2);
	C = [zeros(1,n_in-1) 1 zeros(1,n_im-1) 1];
	% Compute F = C * inv(zI-A) at each z in zSet
	F = zeros(length(zSet),order);
	I = eye(order);
	for i = 1:length(zSet)
	    F(i,:) = C * inv(zSet(i)*I-A);
	end
	B(:,2) = F\L1;
	if rot==1
	    ABCD = [A B(:,2); C 0];
	    for i = 1:order
		phi = angle(ABCD(i,end));
		ABCD(i,:) = ABCD(i,:)*exp(-1i*phi);
		ABCD(:,i) = ABCD(:,i)*exp(1i*phi);
	    end
	    [A B2 C D2] = partitionABCD(ABCD);
	    B(:,2) = B2;
	end
	B(1,1) = abs(B(1,2));
    case 'FF'
	B0 = [1;zeros(order-1,1)];
	B = B0; B(end) = bn;
	% Compute F = inv(zI-A)*B at each z in zSet
	F = zeros(order,length(zSet));
	I = eye(order);
	for i = 1:length(zSet)
	    F(:,i) = inv(zSet(i)*I-A) * B;
	end
	C = L1.'/F;
	if rot==1
	    ABCD = [A B; C 0];
	    for i = 2:order-1
		phi = angle(ABCD(end,i));
		ABCD(i,:) = ABCD(i,:)*exp(1i*phi);
		ABCD(:,i) = ABCD(:,i)*exp(-1i*phi);
	    end
	    [A B C D2] = partitionABCD(ABCD);
	end
	B = [-B0 B];
    case 'PFF'
	B0 = [1;zeros(order-1,1)];
	B = B0; B(n_in+1) = 1;
	% Compute F = inv(zI-A)*B at each z in zSet
	F = zeros(order,length(zSet));
	I = eye(order);
	for i = 1:length(zSet)
	    F(:,i) = inv(zSet(i)*I-A) * B;
	end
	C = L1.'/F;
	if rot==1
	    ABCD = [A B; C 0];
	    for i = 2:order-1
		phi = angle(ABCD(end,i));
		ABCD(i,:) = ABCD(i,:)*exp(1i*phi);
		ABCD(:,i) = ABCD(:,i)*exp(-1i*phi);
	    end
	    [A B C D2] = partitionABCD(ABCD);
	end
	B = [B0 B];
    otherwise
	error(sprintf('%s error. Sorry, form "%s" is not supported.\n', ... 
	    mfilename, form));
end

ABCD = [A B; C D];
