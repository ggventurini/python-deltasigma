function i=outsideConvex(x,n,o,tol) 
% function i=outsideConvex(x,n,o,tol<0>). Determine which points are outside
% the convex object described by n'*x+o<=0. Positive tol reduces the number
% of outside points.
if nargin<4
    tol=0;
end
if( size(x,2)*size(n,2) < 1e5 )
    i = sign(sum(n'*x + o(ones(size(x,2),1),:)' >tol));
else % Matrices are too big. Break into chunks
    chunk = round(1e5/size(n,2));
    i = zeros(1,size(x,2));
    for j=1:chunk:size(x,2)
	j2 = min(j+chunk-1,size(x,2));
	i(j:j2) = sign(sum(n'*x(:,j:j2) + o(ones(j2-j+1,1),:)' >tol));
    end
end

