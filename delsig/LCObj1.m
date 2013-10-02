function [objective,constraints] = LCObj1(x,param,max_radius,dbg)
%function [objective,constraints] = LCObj1(x,param,max_radius,dbg)
% The objective and constraints function for the initial optimization
% process used to put the roots of the denominator of the LCBP NTF inside 
% the unit circle.

H = LCoptparam2tf(x,param);
objective = 1;			% No actual objective
rmax = max(abs(H.p{:}));
constraints = rmax - max_radius;

if dbg
    fprintf(1,'x = [ ');
    fprintf(1, '%.4f ',x);
    fprintf(1,']\n');
    fprintf(1,'rmax = %f\n\n',rmax);
end
