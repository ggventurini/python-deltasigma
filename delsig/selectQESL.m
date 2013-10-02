function sv = selectQESL(v,sy)
%sv = selectQESL(v,sy)	Select the elements of a 2*M-element DAC 
%to minimize the selection error, subject to the constraint that
%the number of elements allocated to [Ip In Qp Qn ]
%is n= [(M+real(v))/2 (M-real(v))/2 (M+imag(v))/2 (M-imag(v))/2].
%As a result of these constraints, sum(sy)=v and M elements are 
%allocated to I and Q at all times.

% See ADI Notebook #17, pages 22 & 24
sy = sy(:);

M = length(sy)/2;
n = zeros(1,4);
n(1) = (M+real(v))/2;	n(3) = (M+imag(v))/2;
n(2) = M-n(1);		n(4) = M-n(3);
if any(n<0)
    error('v argument is out of range');
end

sv = zeros(2*M,1);

%Sort the desired usage by real and imaginary parts
[iList iInd] = sort(real(sy));
[qList qInd] = sort(imag(sy));
%Pad the lists with +/-Inf to allow indices to run off the ends
iList = [-Inf; iList; Inf];	iInd = [0; iInd; 0];
qList = [-Inf; qList; Inf];	qInd = [0; qInd; 0];

% Walk through the lists, top-to-bottom (-I,-Q) and bottom-to-top (+I,+Q),
% allocating elements according to the desired usage, 
% if needed for that purpose.
i1 = 2;	i2 = 2*M+1;
q1 = 2;	q2 = 2*M+1;
[preference mi] = max( [iList(i2) -iList(i1) qList(q2) -qList(q1)] );
while preference>=0
    % Determine which element has highest desired usage
    % Use that element the way it wants to be used, if it is needed there.
    % If an element is used, remove it from the lists.
    % Then move on to the next element on the list.
    i = 0; q = 0;
    switch mi
    case 1	% +I
	if n(1)>0
	    sv(iInd(i2)) = 1;
	    n(1) = n(1)-1;
	    i = i2;
	    q = find(qInd == iInd(i));
	end
	i2 = i2-1;
    case 2	% -I
	if n(2)>0
	    sv(iInd(i1)) = -1;
	    n(2) = n(2)-1;
	    i = i1;
	    q = find(qInd == iInd(i));
	end
	i1 = i1+1;
    case 3	% +Q
	if n(3)>0
	    sv(qInd(q2)) = 1j;
	    n(3) = n(3)-1;
	    q = q2;
	    i = find(iInd == qInd(q));
	end
	q2 = q2-1;
    case 4	% -Q
	if n(4)>0
	    sv(qInd(q1)) = -1j;
	    n(4) = n(4)-1;
	    q = q1;
	    i = find(iInd == qInd(q));
	end
	q1 = q1+1;
    end
    if i	% Delete the chosen element from the lists
	iList = iList([1:i-1 i+1:end]);
	iInd = iInd([1:i-1 i+1:end]);
	qList = qList([1:q-1 q+1:end]);
	qInd = qInd([1:q-1 q+1:end]);
	if i1>i;  i1=i1-1; end
	if i2>=i; i2=i2-1; end
	if q1>q;  q1=q1-1; end
	if q2>=q; q2=q2-1; end
    end
    if q2<q1 & i2<i1
	preference = -1;
    else
	[preference mi] = max( [iList(i2) -iList(i1) qList(q2) -qList(q1)] );
    end
end


% Now that elements have been allocated the way they prefer to be,
% two of the constraints should be satisfied
if n(1)+n(2)>0 & n(3)+n(4)>0	% Two constraints not satisfied
				% A trade-off is needed
    if n(1)
	ni = 1;
	i_val = 1;
	iList = -iList;
    else
	ni = 2;
	i_val = -1;
	iList = iList(end:-1:1);
	iInd = iInd(end:-1:1);
    end
    if n(3)
	nq = 3;
	q_val = 1j;
	qList = -qList;
    else
	nq = 4;
	q_val = -1j;
	qList = qList(end:-1:1);
	qInd = qInd(end:-1:1);
    end
    % Examine elements in the iList and qList
    i1 = 2;	q1 = 2;
    while n(1)+n(2)>0 & n(3)+n(4)>0
	if iList(i1) > qList(q1)
	    % Assign the element to Q 
	    sv(iInd(i1)) = q_val;
	    n(nq) = n(nq) -1;
	    i = i1;
	    q = find(qInd == iInd(i));
	else
	    % Assign the element to I 
	    sv(qInd(i1)) = i_val;
	    n(ni) = n(ni) -1;
	    q = q1;
	    i = find(iInd == qInd(q));
	end
	% Delete the chosen element from the lists
	iList = iList([1:i-1 i+1:end]);
	iInd = iInd([1:i-1 i+1:end]);
	qList = qList([1:q-1 q+1:end]);
	qInd = qInd([1:q-1 q+1:end]);
	if i1>i;  i1=i1-1; end
	if q1>q;  q1=q1-1; end
    end
end

% Three constraints should be satisfied.
% Just assign remaining elements where they are needed
ind = iInd( find(iInd ~=0) );
if n(1)
    sv(ind) = 1;
elseif n(2)
    sv(ind) = -1;
elseif n(3)
    sv(ind) = 1j;
elseif n(4)
    sv(ind) = -1j;
end
return
