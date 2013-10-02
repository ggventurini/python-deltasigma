function H=evalTFP(Hs,Hz,f)
% H=evalTFP(Hs,Hz,f)
% Compute the value of a transfer function product Hs*Hz at a frequency f,
% where Hs is a cts-time TF and Hz is a discrete-time TF.
% Both Hs and Hz are SISO zpk objects.
% This function attempts to cancel poles in Hs with zeros in Hz.

szeros = Hs.z{1};
spoles = Hs.p{1};
zzeros = Hz.z{1};
zpoles = Hz.p{1};

slim = min(1e-3,max(1e-5,eps^(1/(1+length(spoles)))));
zlim = min(1e-3,max(1e-5,eps^(1/(1+length(zzeros)))));

H = zeros(size(f));
w = 2*pi*f;	s = j*w;	z=exp(s);
for i=1:length(f)
    wi = w(i);	si = s(i);	zi = z(i);
    if isempty(spoles)
        cancel = 0;
    else
        cancel = abs(si-spoles)<slim;
    end
    if ~cancel
    	% wi is far from a pole, so just use the product Hs*Hz
    	H(i) = evalTF(Hs,si) * evalTF(Hz,zi);
    else
    	% cancel pole(s) of Hs with corresponding zero(s) of Hz
		cancelz = abs(zi-zzeros)<zlim;
		if sum(cancelz) > sum(cancel)
		    H(i) = 0;
		elseif sum(cancelz) < sum(cancel)
	    	    H(i) = Inf;
		else
		    H(i) =  evalRPoly(szeros,si,Hs.k) * ...
			    zi^sum(cancel) * evalRPoly(zzeros(~cancelz),zi,Hz.k) / ... 
			    (evalRPoly(spoles(~cancel),si,1)*evalRPoly(zpoles,zi,1));
		end
    end
end
