function sv = ds_therm(v,M)
% sv = ds_therm(v,M)	sv is an M by length(v) matrix of +/-1s
% where sum(sv)=v and the -1s are located in the bottom rows
sv = -ones(M,length(v));
for i=1:length(v)
    sv(1:(M+v(i))/2,i) = 1;
end
