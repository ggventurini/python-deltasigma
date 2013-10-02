function zpk2 = cancelPZ(zpk1, tol)
% zpk2 = cancelPZ(zpk1, tol=1e-6) Cancel zeros/poles in a zpk system
if nargin<2
	tol = 1e-6;
end

z = zpk1.z{1};
p = zpk1.p{1};
zpk2 = zpk1;	% Captures k, and whether tf is in s or z
% Need to go in reverse order because z gets smaller with each cancellation
for i = length(z):-1:1
	d = z(i) -p;
	cancel = find(abs(d)<tol);
	if cancel
		p(cancel(1)) = [];
		z(i) =[];
	end
end
zpk2.z = z;
zpk2.p = p;
return
