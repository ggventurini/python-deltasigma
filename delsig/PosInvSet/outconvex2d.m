function out = outconvex2d(x,p)
%function out = outconvex2d(x,p)
%Test if each of the x points are inside the convex polygon p,
%and return the number of inequalities failed by each point.
%p is a 2xn counter-clockwise list of vertices,
%with the first vertex duplicated.
n = size(p,2);

% form A,B such that internal points satisfy Ax <= B
A = [	p(2,2:n) - p(2,1:n-1);
	p(1,1:n-1) - p(1,2:n)]';
B = zeros(n-1,1);
for i=1:n-1
    B(i) = A(i,:)*p(:,i);
end

out = sum( A*x > B(:,ones(1,size(x,2))) );
