function [f1_saved,f2_saved,info]=designHBF(fp,delta,debug)
%function [f1,f2,info]=designHBF(fp=0.2,delta=1e-5,debug=0)
%Design a half-band filter which can be realized without general multipliers.
%The filter is a composition of a prototype and sub- filter.
%Input
% fp	The normalized cutoff frequency of the filter. Due to the
%	symmetry imposed by a HBF, the stopband begins at 0.5-fp.
% delta	The absolute value of the deviation of the frequency response from 
%	the ideal values of 1 in the passband and  0 in the stopband.
%
%Output
% f1,f2	The coefficients of the prototype and sub-filters
%	and their canonical-signed digit (csd) representation.
% info	A vector containing the following data (only set when debug=1):
%	complexity	The number of additions per output sample.
%	n1,n2		The length of the f1 and f2 vectors.
%	sbr		The achieved stob-band attenuation (dB).
%	phi		The scaling factor for the F2 filter.

% To Do: Clean up the code a bit more, esp. wrt the use of the struct. arrays.
%	 Use the phi variable to cut down on the number of adders in F2.
%	 Apply a simulated annealing/genetic optimization alg instead
%	 of the ad hoc one I have now.

%Handle the input arguments
parameters = ['fp   ';'delta';'debug'];
defaults = [ 0.2 1e-5 0];
for i=1:length(defaults)
    if i>nargin
       eval([parameters(i,:) '=defaults(i);'])
    elseif eval(['any(isnan(' parameters(i,:) ')) | isempty(' parameters(i,:) ')']) 
       eval([parameters(i,:) '=defaults(i);'])
    end
end

%Try several different values for the fp1 parameter.
%The best values are usually around .04
%Surrender if 3 successive attempts yield progressively greater complexity.
lowest_complexity = Inf;	prev_complexity = Inf;
for fp1 = [.03 .035 .025 .040 .020 .045 .015 .05]
    failed = 0;
    [f1 zetap phi] = designF1( delta, fp1 );
    if zetap == 1	% designF1 failed
	failed = 1;
	if debug
	    fprintf(2,'designF1 failed at fp1=%f\n',fp1);
	end
    end
    if ~failed
	f2 = designF2( fp, zetap, phi );
	n1 = length(f1);	n2 = length(f2);
	if n2 == 0		% designF2 failed
	    failed = 1;
	    if debug
		fprintf(2,'designF2 failed when zetap=%f, phi=%f\n',zetap,phi);
	    end
	end
    end
    if ~failed
	% complexity(+ performance)  = the number of two-input adders (+ sbr)
	complexity =  size([f1.csd],2) + (2*n1-1)*(n2+size([f2.csd],2)-1);
	if debug
	    msg = sprintf('%d adders: n1=%d, n2=%d, (fp1=%.2f, zetap=%.3f, phi=%4.2f)', ...
		complexity, n1, n2, fp1, zetap, phi );
	else
	    msg = '';
	end
	[fresp pbr sbr] = frespHBF([], f1, f2, phi, fp, msg);
	if pbr <= delta & sbr <= delta          
	    complexity = complexity + sbr;
	    if complexity < prev_complexity
		worse = 0;
		if complexity < lowest_complexity 
		    lowest_complexity = complexity;
		    f1_saved = f1;	f2_saved = f2;
		    phi_saved = phi;
		    if debug
			fprintf( 1, '%s\n', msg )
		    end
		end
	    else
		worse = worse + 1;
		if worse > 2
		    break;
		end
	    end
	    prev_complexity = complexity;
	end	    % if pbr <= delta
    end	    
end	    % for fp1

if isinf(lowest_complexity)
    fprintf(1,'%s: Unable to meet the design requirements.\n', mfilename);
elseif debug 
    complexity = floor(lowest_complexity);
    msg = sprintf( 'Final Design: %d adders', complexity);
    [junk pbr sbr] = frespHBF([], f1_saved, f2_saved, phi_saved, fp, msg);
    n1 = length(f1_saved);	n2 = length(f2_saved);
    fprintf(1,'%s (%d,%d,%.0fdB)\n', msg,n1,n2,dbv(sbr));
    info = [ complexity n1 n2 dbv(sbr) phi_saved ];
end
return


function [f1_saved,zetap,phi] = designF1(delta, fp1)
% [f1 zetap phi] = designF1(delta, fp1)		Design the F1 sub-filter
% of a Saramaki halfband filter. This function is called by designHBF.m.
%
% f1    A structure array containing the F1 filter coefficents and
%       Their CSD representation.
% phi	The scaling factor for the F2 filter (imbedded in the f1 coeffs.)

passband = exp(4*pi*j*linspace(0,fp1));
ok = 0;
for n1 = 1:2:7 	% Odd values only
    if n1 == 1
	h = [0.5 0.5];
    else
	h = firpm(2*n1-1,[0 4*fp1 1 1],[1 1 0 0]);
	if ~(abs(sum(h)-1) < 1e-3 )		% remez bug! Use firls instead
	    h = firls(2*n1-1,[0 4*fp1 1-1e-6 1],[1 1 0 0]);
	end
    end
    fresp = abs( polyval(h,passband) );
    if max( abs(fresp-1) ) <= delta
	ok = 1;
	break
    end
end
if ~ok
    zetap = 1;	% Use this as an indication that the function failed.
    return
end

% Transform h(n) to a chebyshev polynomial f1(n)
% Sum(f1(i)*cos(w)^n)|i=1:n1 + Sum(h(n1+i))*cos(n*w))|i=1:n1, n = 2*i-1;
w = pi*rand(1,n1);
cos_w = cos(w);
A = zeros(n1,length(w));
B = zeros(1,n1);
for i = 1:n1
    n = 2*i-1;
    A(i,:) = cos_w .^ n;
    B = B + h(n1+i)* cos(n*w);
end
f1 = B/A;

% Matlab Ver. 5 change:
phivecb = [];

% Optimize the quantized version of f1 to maximize the stopband width 
% ( = acos(zetap) )

zetap = 1;
testPoints = [0 logspace(-2,0,128)] - 1;
for nsd = 3:8
    f1a = f1'; f1b = f1'; 		% First try the unperturbed filter.
    for phia = 1 ./ [1 f1]
	phia = phia / 2^nextpow2(phia); % keep phi in (0.5,1]
	% Try a bunch of coefficients in the current neighborhood,
	% shrinking the neighborhood once 10 successive trial values show no
	% improvement.  If 2 successive shrinkages do no good, try a higher nsd.
	count = 0;
	nohelp = 0;
	neighborhood = .05;
	while neighborhood > 1e-5
	    phivec = phia .^ [1:2:2*n1-1]';
% Matlab Ver. 5 change:
	    if isempty(phivecb); phivecb = phivec; end
	    f1q = bquantize( f1a.*phivec, nsd );
	    F1 = evalF1( [f1q.val], testPoints, phia );
	    fi = find( abs(F1) > delta ); 
	    zeta = -testPoints( max( fi(1)-1, 1 ) );
	    %fprintf(2,'nsd=%d, nbhd= %f, count=%d, zeta = %f, phia=%f\n', ...
	    %  nsd, neighborhood, count, zeta, phia );
	    if zeta < zetap
		count = 0;
		nohelp = 0;
		zetap = zeta;
		f1b = [f1q.val]';
		f1_saved = f1q;
		phi = phia;
		phivecb = phivec;
	    else
		count = count + 1;
	    end
	    if count > 10
		count = 0;
		neighborhood = neighborhood/2;
		nohelp = nohelp +1;
		if nohelp > 2
		    break;
		end
	    end
	    f1a = f1b./phivecb + neighborhood*(rand(size(f1b))-0.5);
	    phia = phia + neighborhood*(rand(1,1)-0.5);
	end
	if zetap < 1	% Found a filter with adequate attn.
	    break;
	end
    end			% for phia ...
    if zetap < 1	% Found a filter with adequate attn.
	break;
    end
end
return

function f2 = designF2(fp,zetap,phi)
% f2 = designF2(fp,zetap,phi)		Design the F2 sub-filter
% of a Saramaki halfband filter.  This function is called by designHBF.m.

% subfilter design:
%   1 - delta2' < |F2/phi| < 1 	for f in [0 fp];
%  -1 < |F2/phi| < -1 + delta2'	for f in [0.5-fp, 0.5];
%   1-delta2' = (1-delta2)/(1+delta2)

delta2 = (1-zetap)/(1+zetap);
%delta2p = 1 - (1-delta2)/(1+delta2);

% determine the minimum order required by the filter
passband = exp(j*linspace(0,4*pi*fp));
for nsub = 3:2:17
    h2 = firpm(nsub,[0 4*fp 1 1], [1 1 0 0]);
    mag = abs( polyval(h2,passband) );
    if max(abs(mag-1)) < delta2;
	break;
    end
end
n2min = (nsub+1)/2;

% Search all n2,nsd pairs, in order of the product n2*(nsd+1)
% allowing fp to be a variable?
success = 0;
nsdmin = 3;	nsdmax = 6;
for product = (nsdmin+1)*n2min:(nsdmax+1)*n2min
    for nsd = nsdmin:nsdmax
    	n2 = product/(nsd+1);
	if floor(n2) ~= n2	% Only take integer n2,nsd pairs
	    break
	end
	nsub = 2*n2-1;
	% Could try a bunch of fp values
	%fprintf(2,'designF2: Trying (n2,nsd2,fp)=(%2d,%2d,%6.4f)\n',n2,nsd,fp);
	h2 = firpm(nsub,[0 4*fp 1 1], [1 1 0 0]);
	h2 =  h2/(phi*(1+delta2));		% Adjust the coefficients.
	f2 = bquantize( h2(n2+1:nsub+1), nsd );
	h2 = (1+delta2)*phi*[f2(n2:-1:1).val f2.val];
	mag = abs( polyval(h2,passband) );
	if max(abs(mag-1)) < delta2;
	    success =1;
	    break;
	end
    end
    if success
	break;
    end
end

if ~success
    f2 = []; 
    q2 = [];
end
return
