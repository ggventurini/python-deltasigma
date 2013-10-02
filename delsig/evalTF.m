function h = evalTF(tf,z)
%h = evalTF(tf,z)
%Evaluates the rational function described by the struct tf
% at the point(s) given in the z vector.
% TF must be either a zpk object or a struct containing
%   form		'zp' or 'coeff'
%   zeros,poles,k	if form=='zp'
%   num,den		if form=='coeff'
%
% In Matlab 5, the ss/freqresp() function does nearly the same thing.

if isobject(tf)		% zpk object
    if strcmp(class(tf),'zpk')
	h = tf.k * evalRPoly(tf.z{1},z) ./ evalRPoly(tf.p{1},z);
    else
	fprintf(1,'%s: Only zpk objects supported.\n', mfilename);
    end
elseif any(strcmp(fieldnames(tf),'form'))
    if strcmp(tf.form,'zp')
	h = tf.k * evalRPoly(tf.zeros,z) ./ evalRPoly(tf.poles,z);
    elseif strcmp(tf.form,'coeff')
	h = polyval(tf.num,z) ./ polyval(tf.den,z);
    else
	fprintf(1,'%s: Unknown form: %s\n', mfilename, tf.form);
    end
else	% Assume zp form
    h = tf.k * evalRPoly(tf.zeros,z) ./ evalRPoly(tf.poles,z);
end
