function ABCD = stuffABCD(a,g,b,c,form)
%ABCD = stuffABCD(a,g,b,c,form='CRFB')
%Compute the ABCD matrix for the specified structure.
%See realizeNTF.m for a list of supported structures.
%mapABCD is the inverse function.

% Code common to all structures.
if nargin<5
    if nargin<4
	fprintf(1,'%s error: Insufficient arguments (%.0f).\n',mfilename,nargin)
	return
    else
	form = 'CRFB';
    end
end

order = length(a);
odd = rem(order,2); even = ~odd;
ABCD = zeros(order+1,order+2);
if length(b)==1
    b = [b zeros(1,order)];
end

switch form
case 'CRFB'
    %C=(0 0...c_n)
    % This is done as part of the construction of A, below
    %B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
    ABCD(:,order+1) = b';
    %B2 = -(a_1 a_2... a_n)
    ABCD(1:order,order+2) = -a';
    diagonal = 1:(order+2):order*(order+1);
    ABCD(diagonal) = ones(1,order);
    subdiag = diagonal((1+even):2:order)+1;
    ABCD(subdiag)= c(1+even:2:order);
    supdiag = subdiag((1+odd):length(subdiag))-2;
    ABCD(supdiag) = -g;
    dly = (2+odd):2:order;	% row numbers of delaying integrators
    ABCD(dly,:) = ABCD(dly,:) + diag(c(dly-1))*ABCD(dly-1,:);

case 'CRFF'
    %B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
    ABCD(:,order+1) = b';
    %B2 = -(c_1 0... 0)
    ABCD(1,order+2) = -c(1);
    diagonal = 1:(order+2):order*(order+1);	% # of elements = order
    ABCD(diagonal) = ones(1,order);
    subdiag = diagonal(1:2:order-1)+1;
    ABCD(subdiag)= c(2:2:order);
    if even
	multg = 1:2:order;	% rows to have g*(following row) subtracted.
	ABCD(multg,:) = ABCD(multg,:) - diag(g)*ABCD(multg+1,:);
    elseif order >= 3
	supdiag = diagonal(3:2:order)-1;
	ABCD(supdiag) = -g;
    end
    multc = 3:2:order;		% rows to have c*(preceding row) added.
    ABCD(multc,:) = ABCD(multc,:) + diag(c(multc))*ABCD(multc-1,:);
    ABCD(order+1,1:2:order) = a(1:2:order);
    for i = 2:2:order
	ABCD(order+1,:) = ABCD(order+1,:) + a(i)*ABCD(i,:);
    end

case 'CIFB'
    %C=(0 0...c_n)
    % This is done as part of the construction of A, below
    %B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
    ABCD(:,order+1) = b';
    %B2 = -(a_1 a_2... a_n)
    ABCD(1:order,order+2) = -a';
    diagonal = 1:(order+2):order*(order+1);
    ABCD(diagonal) = ones(1,order);
    subdiag = diagonal(1:order)+1;
    ABCD(subdiag)= c;
    supdiag = diagonal((2+odd):2:order)-1;
    ABCD(supdiag) = -g;

case 'CIFF'
    %B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
    ABCD(:,order+1) = b';
    %B2 = -(c_1 0... 0)
    ABCD(1,order+2) = -c(1);
    diagonal = 1:(order+2):order*(order+1);
    ABCD(diagonal) = ones(1,order);
    subdiag = diagonal(1:order-1)+1;
    ABCD(subdiag)= c(2:end);
    %C = (a_1 a_2... a_n)
    ABCD(order+1,1:order) = a(1:order);
    supdiag = diagonal((2+odd):2:order)-1;
    ABCD(supdiag) = -g;

case 'CRFBD'
    %C=(0 0...c_n)
    ABCD(order+1,order) = c(order);
    %B1 = (b_1 b_2... b_n), D=(b_n+1 0)
    ABCD(:,order+1) = b';
    %B2 = -(a_1 a_2... a_n)
    ABCD(1:order,order+2) = -a';
    diagonal = 1:(order+2):order*(order+1);
    ABCD(diagonal) = ones(1,order);
    dly = (1+odd):2:order;	% row numbers of delaying integrators
    subdiag = diagonal(dly)+1;
    ABCD(subdiag)= c(dly);
    supdiag = subdiag((1+odd):length(subdiag))-2;
    ABCD(dly,:) = ABCD(dly,:) - diag(g)*ABCD(dly+1,:);
    if order>2
	coupl = 2+even:2:order;
	ABCD(coupl,:) = ABCD(coupl,:) + diag(c(coupl-1))*ABCD(coupl-1,:);
    end

case 'CRFFD'
    diagonal = 1:(order+2):order*(order+1);
    subdiag = diagonal(1:order-1)+1;
    supdiag = diagonal(2:order)-1;
    %B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
    ABCD(:,order+1) = b';
    %B2 = -(c_1 0... 0)
    ABCD(1,order+2) = -c(1);
    ABCD(diagonal) = ones(1,order);
    multc = 2:2:order;		% rows to have c*(preceding row) added.
    if order>2
	ABCD(subdiag(2:2:end)) = c(3:2:end);
    end
    if even
	ABCD(supdiag(1:2:end)) = -g;
    else
	% subtract g*(following row) from the multc rows
	ABCD(multc,:) = ABCD(multc,:) - diag(g)*ABCD(multc+1,:);
    end
    ABCD(multc,:) = ABCD(multc,:) + diag(c(multc))*ABCD(multc-1,:);
    % C
    ABCD(order+1,2:2:order) = a(2:2:order);
    for i = 1:2:order
	ABCD(order+1,:) = ABCD(order+1,:) + a(i)*ABCD(i,:);
    end
    % The above gives y(n+1); need to add a delay to get y(n).
    % Do this by augmenting the states. Note: this means that
    % the apparent order of the NTF is one higher than it acually is.
    [A B C D] = partitionABCD(ABCD,2);
    A = [A zeros(order,1); C 0];
    B = [B; D];
    C = [zeros(1,order) 1];
    D = [0 0];
    ABCD = [A B; C D];

case 'PFF'
    %B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
    ABCD(:,order+1) = b';
    odd_1 = odd;		% !! Bold assumption !!
    odd_2 = 0;			% !! Bold assumption !!
    gc = g.* c(1+odd_1:2:end);
    theta = acos(1-gc/2);	
    if odd_1
	theta0 = 0;
    else
	theta0 = theta(1);
    end
    order_2 = 2*length(find(abs(theta-theta0)>0.5));
    order_1 = order - order_2;
    %B2 = -(c_1 0...0 c_n 0...0)
    ABCD(1,order+2) = -c(1);
    ABCD(order_1+1,order+2) = -c(order_1+1);
    diagonal = 1:(order+2):order*(order+1);	% # of elements = order
    ABCD(diagonal) = ones(1,order);
    i = [1:2:order_1-1 order-order_2+1:2:order]
    subdiag = diagonal(i)+1;
    ABCD(subdiag)= c(i+1);
    if odd_1
	if order_1 >= 3
	    supdiag = diagonal(3:2:order_1)-1;
	    ABCD(supdiag) = -g(1:(order_1-1)/2);
	end
    else
	multg = 1:2:order_1;	% rows to have g*(following row) subtracted.
	ABCD(multg,:) = ABCD(multg,:) - diag(g(1:order_1/2))*ABCD(multg+1,:);
    end
    if odd_2
	if order_2 >= 3
	    supdiag = diagonal(order_1+2:2:order)-1;
	    ABCD(supdiag) = -g(1:(order_1-1)/2);
	end
    else
	multg = order_1+1:2:order; % rows to have g*(following row) subtracted.
	gg = g((order_1-odd_1)/2+1:end);
	ABCD(multg,:) = ABCD(multg,:) - diag(gg)*ABCD(multg+1,:);
    end
    % Rows to have c*(preceding row) added.
    multc = [3:2:order_1 order_1+3:2:order];		
    ABCD(multc,:) = ABCD(multc,:) + diag(c(multc))*ABCD(multc-1,:);
    % C portion of ABCD
    i = [1:2:order_1 order_1+1:2:order];
    ABCD(order+1,i) = a(i);
    for i = [2:2:order_1 order_1+2:2:order]
	ABCD(order+1,:) = ABCD(order+1,:) + a(i)*ABCD(i,:);
    end

case 'Stratos'
% code copied from case 'CIFF':
    %B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
    ABCD(:,order+1) = b';
    %B2 = -(c_1 0... 0)
    ABCD(1,order+2) = -c(1);
    diagonal = 1:(order+2):order*(order+1);
    ABCD(diagonal) = ones(1,order);
    subdiag = diagonal(1:order-1)+1;
    ABCD(subdiag)= c(2:end);
% code based on case 'CRFF':
    multg = 1+odd:2:order-1;	% rows to have g*(following row) subtracted.
    ABCD(multg,:) = ABCD(multg,:) - diag(g)*ABCD(multg+1,:);
% code copied from case 'CIFF':
    %C = (a_1 a_2... a_n)
    ABCD(order+1,1:order) = a(1:order);

otherwise
    fprintf(1,'%s error: Form %s is not yet supported.\n',mfilename,form)
end
