function y = impL1(arg1,n)
% y=impL1(ntf,n=10) 
% Compute the impulse response from the comparator
% output to the comparator input for the given NTF.
% n is the (optional) number of points (10).
% 
% This function is useful when verifying the realization
% of a NTF with a specified topology.
if nargin<2
    n=10;
end
if isobject(arg1) & strcmp(class(arg1),'zpk')
    z = arg1.z{1};
    p = arg1.p{1};
elseif isstruct(arg1)
    if any(strcmp(fieldnames(arg1),'zeros'))
	z = arg1.zeros;
    else
	error('No zeros field in the NTF.')
    end
    if any(strcmp(fieldnames(arg1),'poles'))
	p = arg1.poles;
    else
	error('No poles field in the NTF.')
    end
end

lf_den = padr(poly(z),length(p)+1);
lf_num = lf_den-poly(p);
if any(imag([lf_num lf_den]))
    % Complex loop filter
    lfr_den = real( conv(lf_den,conj(lf_den)) );
    lfr_num = conv(lf_num,conj(lf_den));
    lf_i = tf( real(lfr_num), lfr_den, 1);
    lf_q = tf( imag(lfr_num), lfr_den, 1);
    y = impulse(lf_i,n) + 1i * impulse(lf_q,n);
else
    y = impulse(tf(lf_num,lf_den,1),n);
end
