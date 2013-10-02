function [a,g,b,c] = realizeNTF(ntf,form,stf)
% [a,g,b,c] = realizeNTF(ntf,form='CRFB',stf=1)
% Convert a noise transfer function into coefficients for the desired structure.
% Supported structures are
%	CRFB	Cascade of resonators, feedback form.
% 	CRFF	Cascade of resonators, feedforward form.
%	CIFB	Cascade of integrators, feedback form. 
%	CIFF	Cascade of integrators, feedforward form.
% 	CRFBD	CRFB with delaying quantizer.
% 	CRFFD	CRFF with delaying quantizer.
%       PFF	Parallel feed-forward.
%	Stratos A CIFF-like structure with non-delaying resonator feedbacks,
%               contributed in 2007 by Jeff Gealow
% See the accompanying documentation for block diagrams of each structure
%
% The order of the NTF zeros must be (real, complex conj. pairs).
% The order of the zeros is used when mapping the NTF onto the chosen topology.
%
% stf is a zpk transfer function

% The basic idea is to equate the loop filter at a set of
% points in the z-plane to L1 = 1-1/ntf at those points.
stderr = 2;
% Handle the input arguments
parameters = {'ntf';'form';'stf'};
defaults = {NaN, 'CRFB', []};
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end

% Code common to all functions
ntf_p = ntf.p{1};
ntf_z = ntf.z{1};
order = length(ntf_p);
order2 = floor(order/2);
odd = order - 2*order2;

a = zeros(1,order);
g = zeros(1,order2);
b = zeros(1,order+1);
c = ones(1,order);

% Choose a set of points in the z-plane at which to try to make L1 = 1-1/H
% I don't know how to do this properly, but the code below seems to work 
% for order <= 11 
N = 200;
min_distance = 0.09;
C = zeros(N,1);
j = 1;
for i=1:N
    z = 1.1*exp(2i*pi*i/N);
    if all( abs(ntf_z-z) > min_distance )
	C(j) = z;
	j = j+1;
    end
end
C(j:end) = [];
zSet = C;

switch form
case 'CRFB'
    %Find g
    %Assume the roots are ordered, real first, then cx conj. pairs 
    for i=1:order2
	g(i) = 2*(1-real(ntf_z(2*i-1+odd)));
    end
    L1 = zeros(1,order);
    %Form the linear matrix equation a*T*=L1
    for i=1:order*2
	z = zSet(i);
	%L1(z) = 1-1/H(z)
	L1(i) = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z);
	Dfactor = (z-1)/z;
	product=1;
	for j=order:-2:(1+odd)
	    product = z/evalRPoly(ntf_z((j-1):j),z)*product;
	    T(j,i) = product*Dfactor;
	    T(j-1,i) = product;
	end
	if( odd )
	    T(1,i)=product/(z-1);
	end
    end
    a = -real(L1/T);
    if isempty(stf)
	b = a;
	b(order+1) = 1;
    end

case 'CRFF'
    %Find g
    %Assume the roots are ordered, real first, then cx conj. pairs
    for i=1:order2
	g(i) = 2*(1-real(ntf_z(2*i-1+odd)));
    end
    L1 = zeros(1,order);
    %Form the linear matrix equation a*T*=L1
    for i=1:order*2
	z = zSet(i);
	%L1(z) = 1-1/H(z)
	L1(i) = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z);
	if( odd )
	    Dfactor = z-1; 
	    product = 1/Dfactor;
	    T(1,i) = product;
	else
	    Dfactor = (z-1)/z;
	    product=1;
	end
	for j=1+odd:2:order
	    product = z/evalRPoly(ntf_z(j:j+1),z)*product;
	    T(j,i) = product*Dfactor;
	    T(j+1,i) = product;
	end
    end
    a = -real(L1/T);
    if isempty(stf)
	b = [ 1 zeros(1,order-1) 1];
    end

case 'CIFB'
    %Assume the roots are ordered, real first, then cx conj. pairs
    %Note ones which are moved significantly.
    if any( abs(real(ntf_z)-1) > 1e-3)
	fprintf(stderr,'%s Warning: The ntf''s zeros have had their real parts set to one.\n', mfilename);
    end
    ntf_z = 1 + 1i*imag(ntf_z);
    for i=1:order2
	g(i) = imag(ntf_z(2*i-1+odd))^2;
    end
    L1 = zeros(1,order);
    %Form the linear matrix equation a*T*=L1
    for i=1:order*2
	z = zSet(i);
	%L1(z) = 1-1/H(z)
	L1(i) = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z);
	Dfactor = (z-1);
	product = 1;
	for j=order:-2:(1+odd)
	    product = product/evalRPoly(ntf_z((j-1):j),z);
	    T(j,i) = product*Dfactor;
	    T(j-1,i) = product;
	end
	if( odd )
	    T(1,i) = product/(z-1);
	end
    end
    a = -real(L1/T);
    if isempty(stf)
	b = a;
	b(order+1) = 1;
    end

case 'CIFF'
    %Assume the roots are ordered, real first, then cx conj. pairs
    %Note ones which are moved significantly.
    if any( abs(real(ntf_z)-1) > 1e-3 )
	fprintf(stderr,'%s Warning: The ntf''s zeros have had their real parts set to one.\n', mfilename);
    end
    ntf_z = 1 + 1i*imag(ntf_z);
    for i=1:order2
	g(i) = imag(ntf_z(2*i-1+odd))^2;
    end
    L1 = zeros(1,order);
    %Form the linear matrix equation a*T*=L1
    for i=1:order*2
	z = zSet(i);
	%L1(z) = 1-1/H(z)
	L1(i) = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z);
	Dfactor = (z-1);
	if( odd )
	    product = 1/(z-1);
	    T(1,i) = product;
	else
	    product = 1;
	end
	for j=odd+1:2:order-1
	    product = product/evalRPoly(ntf_z(j:j+1),z);
	    T(j,i) = product*Dfactor;
	    T(j+1,i) = product;
	end
    end
    a = -real(L1/T);
    if isempty(stf)
	b = [ 1 zeros(1,order-1) 1];
    end

case 'CRFBD'
    %Find g
    %Assume the roots are ordered, real first, then cx conj. pairs
    for i=1:order2
	g(i) = 2*(1-real(ntf_z(2*i-1+odd)));
    end
    L1 = zeros(1,order);
    %Form the linear matrix equation a*T*=L1
    for i=1:order*2
	z = zSet(i);
	%L1(z) = 1-1/H(z)
	L1(i) = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z);
	Dfactor = (z-1);
	product=1/z;
	for j=order:-2:(1+odd)
	    product = z/evalRPoly(ntf_z((j-1):j),z)*product;
	    T(j,i) = product*Dfactor;
	    T(j-1,i) = product;
	end
	if( odd )
	    T(1,i)=product*z/(z-1);
	end
    end
    a = -real(L1/T);
    if isempty(stf)
	b = a;
	b(order+1) = 1;
    end

case 'CRFFD'
    %Find g
    %Assume the roots are ordered, real first, then cx conj. pairs
    for i=1:order2
	g(i) = 2*(1-real(ntf_z(2*i-1+odd)));
    end
    %zL1 = z*(1-1/H(z))
    zL1 = zSet .* (1-1./evalTF(ntf,zSet));
    %Form the linear matrix equation a*T*=zL1
    for i=1:order*2
	z = zSet(i);
	if( odd )
	    Dfactor = (z-1)/z;
	    product = 1/Dfactor;
	    T(1,i) = product;
	else
	    Dfactor = z-1; 
	    product=1;
	end
	for j=1+odd:2:order
	    product = z/evalRPoly(ntf_z(j:j+1),z)*product;
	    T(j,i) = product*Dfactor;
	    T(j+1,i) = product;
	end
    end
    a = -real(zL1/T);
    if isempty(stf)
	b = [ 1 zeros(1,order-1) 1];
    end

case 'PFF'
    %Find g
    %Assume the roots are ordered, real first, then cx conj. pairs
    % with the secondary zeros after the primary zeros
    for i=1:order2
	g(i) = 2*(1-real(ntf_z(2*i-1+odd)));
    end
    % Find the dividing line between the zeros
    theta0 = abs(angle(ntf_z(1)));
    % !! 0.5 radians is an arbitrary separation !!
    i = find( abs(abs(angle(ntf_z)) - theta0) > 0.5 );
    order_1 = i(1)-1;
    order_2 = order-order_1;
    if length(i) ~= order_2
	keyboard
	error('For the PFF form, the NTF zeros must be sorted into primary and secondary zeros');
    end
    odd_1 = mod(order_1,2);
    odd_2 = mod(order_2,2);
    L1 = zeros(1,order);
    %Form the linear matrix equation a*T*=L1
    for i=1:order*2
	z = zSet(i);
	%L1(z) = 1-1/H(z)
	L1(i) = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z);
	if( odd_1 )
	    Dfactor = z-1; 
	    product = 1/Dfactor;
	    T(1,i) = product;
	else
	    Dfactor = (z-1)/z;
	    product=1;
	end
	for j=1+odd_1:2:order_1
	    product = z/evalRPoly(ntf_z(j:j+1),z)*product;
	    T(j,i) = product*Dfactor;
	    T(j+1,i) = product;
	end
	if( odd_2 )
	    Dfactor = z-1; 
	    product = 1/Dfactor;
	    T(order_1+1,i) = product;
	else
	    Dfactor = (z-1)/z;
	    product=1;
	end
	for j=order_1+1+odd_2:2:order
	    product = z/evalRPoly(ntf_z(j:j+1),z)*product;
	    T(j,i) = product*Dfactor;
	    T(j+1,i) = product;
	end
    end
    a = -real(L1/T);
    if isempty(stf)
	b = [ 1 zeros(1,order_1-1) 1 zeros(1,order_2-1) 1];
    end

case 'Stratos'
% code copied from case 'CRFF':
    %Find g
    %Assume the roots are ordered, real first, then cx conj. pairs
    for i=1:order2
	g(i) = 2*(1-real(ntf_z(2*i-1+odd)));
    end
% code copied from case 'CIFF':
    L1 = zeros(1,order);
    %Form the linear matrix equation a*T*=L1
    for i=1:order*2
	z = zSet(i);
	%L1(z) = 1-1/H(z)
	L1(i) = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z);
	Dfactor = (z-1);
	if( odd )
	    product = 1/(z-1);
	    T(1,i) = product;
	else
	    product = 1;
	end
	for j=odd+1:2:order-1
	    product = product/evalRPoly(ntf_z(j:j+1),z);
	    T(j,i) = product*Dfactor;
	    T(j+1,i) = product;
	end
    end
    a = -real(L1/T);
    if isempty(stf)
	b = [ 1 zeros(1,order-1) 1];
    end

end

if ~isempty(stf)
    % Compute the TF from each feed-in to the output 
    % and solve for coefficients which yield the best match
    % THIS CODE IS NOT OPTIMAL, in terms of computational efficiency.
    stfList = cell(1,order+1);
    for i = 1:order+1
		bi = zeros(1,order+1); bi(i)=1;
		ABCD = stuffABCD(a,g,bi,c,form);
		if strcmp(form(3:4),'FF')
		    ABCD(1,order+2) = -1;	% b1 is used to set the feedback
		end
		[junk stfList{i}] = calculateTF(ABCD);
    end
    % Build the matrix equation b A = x and solve it.
    A = zeros(order+1,length(zSet));
    for i = 1:order+1
		A(i,:) = evalTF(stfList{i},zSet);
    end
    x = evalTF(stf,zSet);
    x = x(:).';
    b = real(x/A);
end
