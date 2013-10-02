function y = evalRPoly(roots,x,k)
%function y = evalRPoly(roots,x,k=1)
%Compute the value of a polynomial which is given in terms of its roots.
if(nargin<3)
    k=1;
end

y = k(ones(size(x)));
roots = roots(~isinf(roots));        % remove roots at infinity
for(i=1:length(roots))
    y = y.*(x-roots(i));
end
