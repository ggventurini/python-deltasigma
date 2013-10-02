function [pplus, eplus, pminus, eminus] = dssplit2d(u,ABCD,p)
%function [pplus, eplus, pminus, eminus] = dssplit2d(u,ABCD,p)
%Split a convex polygon p into "plus" and "minus" polygons.
% C * pplus + D1*u >=0, C * pminus + D1*u <=0. 
%p is given as a sequential list of vertices 
%with the first vertex replicated at the end of the list.
%ABCD describes the modulator structure, 
%and u is the modulator input.
%Limitation: D1 must be zero if u is a range.

n = size(ABCD,1)-1;
C = ABCD(n+1, 1:n);
D1= ABCD(n+1, n+1);	% D2=ABCD(n+1,n+2) must be zero
N = size(p,2);
if length(u)==1
    D1u = D1*u;
    y = C*p + D1u(ones(1,N));
else
    if D1 ~= 0
	fprintf('%s: Error. D1 must be zero when u is a range.\n');
	return;
    else
	y = C*p;
    end
end

sign1 = sgn(y(1));
i = find( sgn(y) ~= sign1 );
i1 = i(1);		% First change of sign.
pa = dscut( p(:,i1-1),y(i1-1), p(:,i1),y(i1) );
i2 = i(length(i));	% Second change of sign.
pb = dscut( p(:,i2),y(i2), p(:,i2+1),y(i2+1) );
if sign1 > 0
    pminus = [pa p(:,i) pb pa ];
    pplus = [p(:,1:i1-1) pa pb p(:,i2+1:N)];
else
    pplus = [pa p(:,i) pb pa ];
    pminus = [p(:,1:i1-1) pa pb p(:,i2+1:N)];
end
ne = size(pplus,2);
eplus = [1:ne; [2:ne] 1];
ne = size(pminus,2);
eminus = [1:ne; [2:ne] 1];
