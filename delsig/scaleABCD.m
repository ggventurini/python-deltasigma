function [ABCDs,umax,S]=scaleABCD(ABCD,nlev,f,xlim,ymax,umax,N_sim,N0)
%[ABCDs,umax,S]=scaleABCD(ABCD,nlev=2,f=0,xlim=1,ymax=nlev+5,umax,N=1e5,N0=10)
%Scale the loop filter of a general delta-sigma modulator for dynamic range.
%
%ABCD	The state-space description of the loop filter.
%nlev	The number of levels in the quantizer.
%xlim	A vector or scalar specifying the limit for each state variable.
%ymax	The stability threshold. Inputs that yield quantizer inputs above ymax
%       are considered to be beyond the stable range of the modulator.
%umax	The maximum allowable input amplitude. umax is calculated if it
%	is not supplied.
%
%ABCDs	The state-space description of the scaled loop filter.
%S	The diagonal scaling matrix, ABCDs = [S*A*Sinv S*B; C*Sinv D];
%	xs = S*x;

if nargin<1
    fprintf(1,'Error. scaleABCD needs at least one arguments.\n');
    return
end

% Handle the input arguments
parameters = {'ABCD' 'nlev' 'f' 'xlim' 'ymax' 'umax' 'N_sim' 'N0'};
defaults = { NaN 2 0 1 NaN NaN 1e5 10};
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end
if isnan(ymax)
   ymax=nlev+5;
end
order = size(ABCD,1)-1;
if length(xlim)==1
    xlim = xlim(ones(1,order));	% Convert scalar xlim to a vector
end
if isreal(ABCD)
    quadrature = 0;
else
    quadrature = 1;
end

randn('seed',0);	% So that this function is repeatable

% Envelope for smooth start-up 
raised_cosine = 0.5*(1-cos(pi/N0*[0:N0-1]));
if isnan(umax)
    %Simulate the modulator with DC or sine wave inputs to detect its stable
    %input range. 
    %First get a rough estimate of umax.
    ulist =(0.1:0.1:1.0) * (nlev-1);
    umax = nlev-1;
    N = 1e3;
    u0 = [ exp(2i*pi*f*[-N0:-1]) .* raised_cosine ...
           exp(2i*pi*f*[0:N-1]) ] + 0.01*[1 1i]*randn(2,N+N0);
    if ~quadrature
    	u0 = real(u0);
    end
    for u = ulist
	if ~quadrature
	    [v x xmax y] = simulateDSM(u*u0,ABCD,nlev);
	else
	    [v x xmax y] = simulateQDSM(u*u0,ABCD,nlev);
	end
	if max(abs(y))>ymax
	    umax = u; %umax is the smallest input found which causes 'instability'
	    break
	end
    end
    if umax == ulist(1)
	msg = sprintf('Modulator is unstable even with an input amplitude of %.1f.',umax);
	error(msg);
    end
end
    
%More detailed simulation
N = N_sim;
u0 = [ exp(2i*pi*f*[-N0:-1]) .* raised_cosine ...
       exp(2i*pi*f*[0:N-1]) ] + 0.01*[1 1i]*randn(2,N+N0);
if ~quadrature
	u0 = real(u0);
end
maxima = zeros(1,order)-1;
ulist = linspace(0.7*umax,umax,10);
for u = ulist
	if ~quadrature
		[v x xmax y] = simulateDSM(u*u0,ABCD,nlev);
	else
		[v x xmax y] = simulateQDSM(u*u0,ABCD,nlev);
	end
    if max(abs(y))>ymax
       break
    end
    % We need to do this at least once.
    umax = u;	% umax is the largest input which keeps |y| < ymax
    maxima = max([maxima; xmax']);
end

%Scale the modulator so that all states are at most xlim.
scale = maxima ./ xlim;
S = diag(1./scale); Sinv = diag(scale); % xs = S * x;
[A B C D]= partitionABCD(ABCD);
ABCDs = [S*A*Sinv S*B; C*Sinv D];
