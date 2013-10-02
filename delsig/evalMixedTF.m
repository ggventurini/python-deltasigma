function H=evalMixedTF(tf,f,df)
% H=evalMixedTF(tf,f)
% Compute the mixed transfer function tf at a frequency f.
% tf is a struct array with fields Hs and Hz wich represent 
% a continuous-time/discrete-time tfs which must be multiplied together
% and then added up.
if nargin<3
    df = 1e-5;
end
H = zeros(size(f));
for i = 1:length(tf)
    H = H + evalTFP(tf(i).Hs,tf(i).Hz,f);
end

err = find(isnan(H)|isinf(H));
if ~isempty(err)
    % Need to fill in the holes. !! Use a "nearby" frequency point
    for i=err
	H(i) = evalMixedTF(tf,f(i)+df,df*10);
    end
end
