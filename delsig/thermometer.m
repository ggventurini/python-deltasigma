function t = thermometer(x,m)
% t = thermometer(x,m)	t is an m by length(x) matrix wherein the first
% x(i) components of column i are one;
t = zeros(m,length(x));
for i=1:length(x)
    t(1:x(i),i) = ones(x(i),1);
end
