function sv = ESLselect(v,sy,dw,df)
%sv = ESLselect(v,sy,dw,df)	Select the elements of a multi-element
%DAC to minimize the selection error, subject to the constraint that 
%the nominal DAC output is v, i.e. v = sv'*dw.
%df is a vector of dac error factors. de is defined s.t. df'*de = 0.
%Assume that the preferred usage order is that given by dw.

% Go through sv possibilities one by one, until one which meets the
% v = sv'*dw constraint is found.

if v<0 | v>sum(dw)
    error('v argument is invalid (too large or too small)');
end

n = length(dw);
sv = zeros(n,1);
if v==0; return; end;

[junk possibilities] = sort(-sy); 	%Determine the element priority

% Code to speed up selection in the usual case where all
% DAC elements are weighted at one. Suggested by J. A. Cherry, 5-Mar-97
if dw==ones(n,1);
    sv(possibilities(1:v)) = ones(1,v);
    return;
end;

i = 1;				% Selection level
pointer = 1;			% Array of pointers to selected elements
selected = [];			% Selected elements
while 1
    while pointer(i)>n;
	% backtrack
	i = i-1;
	if i==0
	    break;	%failure!
	end
	pointer(i) = pointer(i)+1;
	selected = selected(1:i);
    end
    selected(i) = possibilities(pointer(i));
    dv = sum(dw(selected));
    if dv==v
	break;		%success!
    elseif dv<v
	% Proceed to the next level of selection
	pointer(i+1) = pointer(i)+1;
	i = i+1;
    else
	% Try the next element at the current level.
	pointer(i) = pointer(i)+1;
    end
end

sv(selected) = ones(1,length(selected));
