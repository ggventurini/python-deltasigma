function dotplot(p,fmt)
%function dotplot(p,fmt)
%Plot dots with format fmt
if nargin < 2 
    fmt = '.';
end
if size(p,1) == 2
    plot( p(1,:), p(2,:), fmt);
elseif size(p,1) == 3
    subplot(221);
    plot(p(1,:),p(2,:),fmt);
    subplot(222);
    plot3(p(1,:),p(2,:),p(3,:),fmt);
    subplot(224);
    plot(p(1,:),p(3,:),fmt);
    subplot(223);
    plot(p(2,:),p(3,:),fmt);
elseif size(p,1) == 4
    subplot(221);
    plot3(p(1,:),p(2,:),p(3,:),fmt);
    subplot(222);
    plot3(p(1,:),p(2,:),p(4,:),fmt);
    subplot(223);
    plot3(p(1,:),p(3,:),p(4,:),fmt);
    subplot(224);
    plot3(p(2,:),p(3,:),p(4,:),fmt);
end
