function [a,g,b,c] = mapABCD(ABCD,form)
%[a,g,b,c] = mapABCD(ABCD,form='CRFB')
%Compute the coefficients for the specified structure.
%See realizeNTF.m for a list of supported structures.
%stuffABCD is the inverse function.

% Code common to all structures.
if nargin<2
    if nargin<1
	fprintf(1,'%s error: Insufficient arguments (%.0f).\n',mfilename,nargin)
	return
    else
	form = 'CRFB';
    end
end

order=size(ABCD,1)-1;
odd = rem(order,2); even = ~odd;
diagonal = 1:(order+2):order*(order+1);
subdiag = diagonal(1:order)+1;
supdiag = diagonal(2+odd:2:order)-1;

%Do the mapping.
switch form
case {'CRFB', 'CIFB', 'CRFBD'}
    c = ABCD(subdiag);
    g = -ABCD(supdiag);
    if strcmp(form,'CRFB') 
	dly = (2+odd):2:order;		%row numbers of delaying integrators.
	ABCD(dly,:) = ABCD(dly,:)-diag(c(dly-1))*ABCD(dly-1,:);
    elseif strcmp(form,'CRFBD') 
	dly = (1+odd):2:order;		%row numbers of delaying integrators.
	ABCD(dly,:) = ABCD(dly,:)+diag(g)*ABCD(dly+1,:);
	if order>2
	    coupl = 2+even:2:order;
	    ABCD(coupl,:) = ABCD(coupl,:) - diag(c(coupl-1))*ABCD(coupl-1,:);
	end
    end
    a = -ABCD(1:order,order+2)';
    b = ABCD(:,order+1)';

case 'CRFF'
    c = [-ABCD(1,order+2) ABCD(subdiag(1:end-1))];
    g = -ABCD(supdiag);
    if even
        multg = 1:2:order;	%Rows to have g*(following row) added.
        ABCD(multg,:) = ABCD(multg,:)+diag(g)*ABCD(multg+1,:);
    end
    multc=3:2:order;	%Rows to have c*(preceding row) subtracted.
    ABCD(multc,:)=ABCD(multc,:)-diag(c(multc))*ABCD(multc-1,:);
    a(2:2:order)=ABCD(order+1,2:2:order);
    for i=2:2:order                   %Recover a coeff.
        ABCD(order+1,:)=ABCD(order+1,:)-a(i)*ABCD(i,:);
    end
    a(1:2:order)=ABCD(order+1,1:2:order);
    b = ABCD(:,order+1)';

case 'CRFFD'
    % CRFFD has an extra order. Correct the common variables
    order = order-1;
    odd = rem(order,2); even = ~odd;
    diagonal = diagonal(1:order);
    subdiag = diagonal(1:order-1)+1;
    supdiag = diagonal(2+odd:2:order)-1;
    g = -ABCD(supdiag);
    c = [-ABCD(1,order+3) ABCD(subdiag)];
    a = zeros(1,order);
    for i=1:2:order                   %Recover the a coefficients
		a(i) = ABCD(order+1,i);
		ABCD(order+1,:) = ABCD(order+1,:) - a(i)*ABCD(i,:);
    end
    a(2:2:order) = ABCD(order+1,2:2:order);
    b = ABCD(1:order+1,order+2)';
    for i=2:2:order                   %Recover the b coefficients
		b(i) = b(i) - c(i) *b(i-1);
		if odd
		    b(i) = b(i) + g(i/2) *b(i+1);
		end
    end
    yscale =  ABCD(order+2,order+1);
    a = a*yscale;
    b(end) = b(end)*yscale;

case {'CIFF','Stratos'}
    a = ABCD(order+1,1:order);
    c = [-ABCD(1,order+2) ABCD(subdiag(1:end-1))];
    g = -ABCD(supdiag);
    b = ABCD(:,order+1)';

otherwise
    fprintf(1,'%s error: Form %s is not yet supported.\n',mfilename,form)
end
