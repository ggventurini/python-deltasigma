function h = polyplot(p,fmt)
%function polyplot(p,fmt)
%Plot a polygon given by the point list p, with format fmt
if nargin < 2 
    fmt = '-';
end
n = size(p,2);
if p(:,n) == p(:,1)
    h = plot( p(1,:), p(2,:), fmt );
else
    h = plot( [p(1,:) p(1,1)], [p(2,:) p(2,1)], fmt );
end
