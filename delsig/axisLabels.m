function s = axisLabels(range,incr)
%function s = axisLabels(range,incr)
range(abs(range)<1e-6) = 0;
s = cell(1,length(range));
if length(incr) == 2
    first = incr(2);
    incr = incr(1);
else
    first = 1;
end
for i=first:incr:length(range)
    s{i}=sprintf('%g',range(i));
end
