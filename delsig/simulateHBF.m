function y = simulateHBF(x,f1,f2,mode)
% y = simulateHBF(x,f1,f2,mode=0)	Simulate a Saramaki half-band filter.
% The f1 and f2 vectors contain coefficients for the structure.
% (f1 and f2 can also be struct arrays like those returned from designHBF.m
% i.e. struct arrays whose .val fields contain the coeffiecients.)
% The mode flag determines whether the input is filtered, interpolated,
% or decimated according to the following table:
%
% mode = 0:	Plain filtering, no interpolation or decimation.
% mode = 1:	The input is interpolated.
% mode = 2:	The output is decimated, even samples are taken.
% mode = 3:	The output is decimated, odd samples are taken.

% Handle the input arguments
if nargin < 1
    fprintf(2, '%s error. Insufficient arguments.\n', mfilename);
    return
end
parameters = {'x' 'f1' 'f2' 'mode'};
defaults = { NaN NaN NaN 0 };
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end
if isnumeric(f1) & isnan(f1)
    [f1 f2] = exampleHBF(4);
end
if isstruct(f1)	% Presume that f1 is a {.val,.csd} struct
    f1 = [f1.val];
end
if isstruct(f2)	% Presume that f2 is a {.val,.csd} struct
    f2 = [f2.val];
end
x = x(:);
f1 = f1(:);		f2 = f2(:);
n1 = length(f1);	n2 = length(f2);
f2imp = [f2(end:-1:1)' f2'];
if mode==0	% Insert zeros
    f2imp = [f2imp; zeros(1,2*n2)];
    f2imp = f2imp(1:4*n2-1);
end
F2 = tf(f2imp,[1 zeros(1,length(f2imp)-1)],1);

switch mode

  case 0	% Plain
    up = lsim(F2,x);
    y = 0.5*delay(x,2*n2-1) + f1(1)*up;
    for i=2:n1
	up = lsim(F2,up);
	up = lsim(F2,up);
	y = f1(i)*up + delay(y,4*n2-2);
    end
    
  case 1	% Interpolating
    up = zeros(size(x));
    nz = 2*n2-1;
    for i=n1:-1:1
	if i==1
	    up = lsim(F2,up+f1(i)*x);
	    nz = n2-1;
	else
	    up = lsim(F2,up+f1(i)*x);
	    up = lsim(F2,up);
	end
	x = delay(x,nz);
    end
    y = [2*up'; x']; y = y(:);	% Interleave the upper and lower streams

  case {2, 3}	% Decimating
    x = x(1:2*floor(end/2));
    if mode==3
	y = 0.5*x(1:2:end);
	up = x(2:2:end);
	nz = n2-1;
    else
	y = 0.5*x(2:2:end);
	up = x(1:2:end);
	nz = n2;
    end
    for i=1:n1
	if i==1
	    up = lsim(F2,up);
	else
	    up = lsim(F2,up);
	    up = lsim(F2,up);
	end
	y = f1(i)*up + delay(y,nz);
	nz = 2*n2-1;
    end

  otherwise
    fprintf(1,'%s: Error. %d is not a valid value for the mode variable.\n', ...
      mfilename, mode );
    return
end
