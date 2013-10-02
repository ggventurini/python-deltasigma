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
