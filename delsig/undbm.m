function v=undbm(p,z) 
% v=undbm(p,z=50) = sqrt(z*10^(p/10-3)) rms voltage equivalent to a power p indBm
if nargin<2
	z = 50;
end
v = sqrt(z*10.^(p/10-3));
