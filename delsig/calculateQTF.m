function [ntf,stf,intf,istf] = calculateQTF(ABCDr)
% [ntf stf intf istf] = calculateQTF(ABCDr) 
% Calculate the noise and signal transfer functions for a quadrature modulator

[A B C D] = partitionABCD(ABCDr,4);
%Construct an ABCD description of the closed-loop system
sys = ss(A+B(:,3:4)*C, B, C, [D(:,1:2) eye(2)], 1);
%Calculate the 2x4 matrix of transfer functions
tfs = tf(sys);
%Reduce these to NTF, STF, INTF and ISTF
if any( tfs(1,3).den{1} ~= tfs(2,4).den{1} )
	error('TF Denominator mismatch. Location 1');
end
ntf_x = tf(0.5*(tfs(1,3).num{1} + tfs(2,4).num{1}), tfs(1,3).den{1}, 1);
intf_x = tf(0.5*(tfs(1,3).num{1} - tfs(2,4).num{1}), tfs(1,3).den{1}, 1);
if any( tfs(1,4).den{1} ~= tfs(2,3).den{1} )
	error('TF Denominator mismatch. Location 2');
end
ntf_y = tf(0.5*(tfs(2,3).num{1} - tfs(1,4).num{1}), tfs(2,3).den{1}, 1);
intf_y = tf(0.5*(tfs(2,3).num{1} + tfs(1,4).num{1}), tfs(2,3).den{1}, 1);
if any( ntf_x.den{1} ~= ntf_y.den{1} )
	error('TF Denominator mismatch. Location 3');
end
if any( tfs(1,1).den{1} ~= tfs(2,2).den{1} )
	error('TF Denominator mismatch. Location 4');
end
stf_x = tf(0.5*(tfs(1,1).num{1} + tfs(2,2).num{1}), tfs(1,1).den{1}, 1);
istf_x = tf(0.5*(tfs(1,1).num{1} - tfs(2,2).num{1}), tfs(1,1).den{1}, 1);
if any( tfs(1,2).den{1} ~= tfs(2,1).den{1} )
	error('TF Denominator mismatch. Location 5');
end
stf_y = tf(0.5*(tfs(2,1).num{1} - tfs(1,2).num{1}), tfs(2,1).den{1}, 1);
istf_y = tf(0.5*(tfs(2,1).num{1} + tfs(1,2).num{1}), tfs(2,1).den{1}, 1);
if any( stf_x.den{1} ~= stf_y.den{1} )
	error('TF Denominator mismatch. Location 6');
end

warning('off'); % suppress warnings about complex TFs
ntf = cancelPZ( zpk( tf(ntf_x.num{1} + 1i* ntf_y.num{1}, ntf_x.den{1}, 1) ) );
intf = cancelPZ( zpk( tf(intf_x.num{1} + 1i* intf_y.num{1}, intf_x.den{1}, 1) ) );
stf = cancelPZ( zpk( tf(stf_x.num{1} + 1i* stf_y.num{1}, ntf_x.den{1}, 1) ) );
istf = cancelPZ( zpk( tf(istf_x.num{1} + 1i* istf_y.num{1}, intf_x.den{1}, 1) ) );
warning('on');

return
