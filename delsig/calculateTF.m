function [ntf,stf] = calculateTF(ABCD,k)
% [ntf,stf] = calculateTF(ABCD,k=1) 
% Calculate the NTF and STF of a delta-sigma modulator whose loop filter
% is described by the ABCD matrix, assuming a quantizer gain of k.
% The NTF and STF are zpk objects.
if nargin < 2 | isnan(k)
    k = 1;
end

[A,B,C,D] = partitionABCD(ABCD);
if size(B,2)>1
    B1 = B(:,1);
    B2 = B(:,2);
else
    B1 = B;
    B2 = B;
end
% Find the noise transfer function by forming the closed-loop
% system (sys_cl) in state-space form.
Acl = A + k*B2*C;
Bcl = [B1 + k*B2*D(1), B2];
Ccl = k*C;
Dcl = [k*D(1) 1];
vn = sscanf(version,'%d');
if vn>6 | all(imag(ABCD)==0)	% real modulator or recent version of MATLAB
    sys_cl = ss(Acl,Bcl,Ccl,Dcl,1);
    tol = min(1e-3,max(1e-6,eps^(1/(size(ABCD,1)))));
    tfs = zpk(sys_cl);
    mtfs = minreal(tfs,tol);
    stf = mtfs(1);
    ntf = mtfs(2);
else	% quadrature modulator and old version of MATLAB
    p = eig(Acl);
    ntfz = eig(A);
    ntf = setPolesAndZeros(ntfz,p,1);
    if nargout>1
        fprintf(1,'Sorry. calculateTF cannot compute the STF for a complex modulator with this version of matlab.\n');
    end
    stf=[];
end
return

% Use a loophole to set complex poles and zeros in zpk objects 
function ztf = setPolesAndZeros(z,p,k)
    tol = 3e-5;		% tolerance for pole-zero cancellation
    Z = zpk(0,[],1,1);
    ztf = zpk([],[],k,1);
    for i = 1:length(p)
	match = abs(p(i)-z)<tol;
	if any(match)
	    f = find(match);
	    z(f(1)) = [];
	else
	    ztf = ztf/(Z-p(i));
	end
    end
    for i = 1:length(z)
	ztf = ztf*(Z-z(i));
    end
return
