function [ABCDc,tdac2] = realizeNTF_ct( ntf, form, tdac, ordering, bp, ABCDc)
% [ABCDc,tdac2] = realizeNTF_ct( ntf, form='FB', tdac, ordering=[1:n], bp=zeros(...), ABCDc)	
% Realize an NTF with a continuous-time loop filter.
%
% Output
% ABCDc      A state-space description of the CT loop filter
%
% tdac2	     A matrix with the DAC timings, including ones
%	     that were automatically added.
%
%
% Input Arguments
% ntf	A noise transfer function in pole-zero form.
%
% form = {'FB','FF'}	
%       A string specifying the topology of the loop filter.
%	For the FB structure, the elements of Bc are calculated
% 	so that the sampled pulse response matches the L1 impulse
% 	respnse.  For the FF structure, Cc is calculated.
%
% tdac	The timing for the feedback DAC(s). If tdac(1)>=1,
% 	direct feedback terms are added to the quantizer.
%	Multiple timings (1 or more per integrator) for the FB 
%       topology can be specified by making tdac a cell array,
%       e.g. tdac = { [1,2]; [1 2]; {[0.5 1],[1 1.5]}; []}; 
%	In this example, the first two integrators have
%	dacs with [1,2] timing, the third has a pair of
% 	dacs, one with [0.5 1] timing and the other with
%	[1 1.5] timing, and there is no direct feedback
% 	DAC to the quantizer
%
% ordering
%	A vector specifying which NTF zero-pair to use in each resonator
% 	Default is for the zero-pairs to be used in the order specified in
%       the NTF.
%
% bp	A vector specifying which resonator sections are bandpass.
%	The default (zeros(...)) is for all sections to be lowpass.
%
% ABCDc The loop filter structure, in state-space form.
%	If this argument is omitted, ABCDc is constructed according 
%       to "form."
% 

% Handle the input arguments
parameters = {'ntf';'form';'tdac';'ordering';'bp';'ABCDc'};
defaults = {NaN, 'FB', [0 1], [], [], []};
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end
ntf_p = ntf.p{1};
ntf_z = ntf.z{1};

order = length(ntf_p);
order2 = floor(order/2);
odd = order - 2*order2;

% compensate for limited accuracy of zero calculation
ntf_z(find(abs(ntf_z - 1) < eps^(1/(1+order)))) = 1;

if iscell(tdac)
    if size(tdac) ~= [order+1 1]
	msg = sprintf(['%s error. For cell array tdac, size(tdac) ' ...
	  'must be [order+1 1].\n'],  mfilename);
	error(msg);
    end
    if form ~= 'FB'
	msg = sprintf(['%s error. Currently only supporting form=''FB'' ' ...
	  'for cell-array tdac'], mfilename); 
	error(msg);
    end
else
    if size(tdac) ~= [1 2]
	msg = sprintf(['%s error. For non cell array tdac, size(tdac) ' ...
	  'must be [1 2].\n'],  mfilename);
	error(msg);
    end
end
if isempty(ordering)
    ordering = [1:order2];
end
if isempty(bp)
    bp = zeros(1,order2);
end
if ~iscell(tdac)
    % Need direct terms for every interval of memory in the DAC
    n_direct = ceil(tdac(2))-1;
    if ((tdac(1)>0) && (tdac(1)<1) && (tdac(2)>1) && (tdac(2)<2))
	n_extra = n_direct-1;     % tdac pulse spans a sample point
    else
	n_extra = n_direct;
    end
    tdac2 = [ -1 -1; 
	       tdac; 
	       0.5*ones(n_extra,1)*[-1 1] + cumsum(ones(n_extra,2),1) ...
                                          + (n_direct - n_extra) ];
else
    n_direct = 0;
    n_extra = 0;
end

if isempty(ABCDc)
    ABCDc = zeros(order+1,order+2);
    % Stuff the A portion
    if odd
	ABCDc(1,1) = real( log( ntf_z(1) ) );
	ABCDc(2,1) = 1;
    end
    for i = 1:order2
	n = bp(i);
	i1 = 2*i + odd - 1;
	zi = 2*ordering(i) + odd - 1;
	w = abs( angle( ntf_z(zi) ) );
	ABCDc(i1+[0 1 2],i1+[0 1]) =[ 0  -w^2
				      1   0
				      n  1-n ];
    end
    ABCDc(1,order+1) = 1;
    ABCDc(1,order+2) = -1;	% 2006.10.02 Changed to -1 to make FF STF have +ve gain at DC
end
Ac = ABCDc(1:order,1:order);
switch form
    case 'FB'
	Cc = ABCDc(order+1,1:order);
	if ~iscell(tdac)
	    Bc = [eye(order) zeros(order,1)];
	    Dc = [zeros(1,order) 1];
	    tp = repmat(tdac,order+1,1);
	else	% Assemble tdac2, Bc and Dc
	    tdac2 = [-1 -1];
	    Bc = [];
	    Dc = [];
	    Bci = [eye(order) zeros(order,1)];
	    Dci = [zeros(1,order) 1];
	    for i=1:length(tdac)
		tdi = tdac{i};
		if iscell(tdi)
		    for j=1:length(tdi)
			tdj = tdi{j};
			tdac2 = [tdac2; tdj];
			Bc = [Bc Bci(:,i)];
			Dc = [Dc Dci(:,i)];
		    end
		elseif ~isempty(tdi)
		    tdac2 = [tdac2; tdi];
		    Bc = [Bc Bci(:,i)];
		    Dc = [Dc Dci(:,i)];
		end
	    end
	    tp = tdac2(2:end,:);
	end
    case 'FF'
	Cc = [eye(order); zeros(1,order)];
   	Bc = [-1; zeros(order-1,1)]; 
	Dc = [zeros(order,1); 1];
	tp = tdac;	% 2008-03-24 fix from Ayman Shabra
    otherwise
	error(sprintf('%s error. Sorry, no code for form "%s".\n', ... 
	    mfilename, form));
end
% Sample the L1 impulse response
n_imp = ceil( 2*order + max(tdac2(:,2)) +1 );
y = impL1(ntf,n_imp);

sys_c = ss( Ac, Bc, Cc, Dc );
yy = pulse(sys_c,tp,1,n_imp,1);
yy = squeeze(yy);
% Endow yy with n_extra extra impulses.
% These will need to be implemented with n_extra extra DACs.
% !! Note: if t1=int, matlab says pulse(sys) @t1 ~=0
% !! This code corrects this problem.
if n_extra>0
    y_right = padb([zeros(1,n_direct); eye(n_direct)], n_imp+1);
    % Replace the last column in yy with an ordered set of impulses
    if (n_direct > n_extra)
	yy = [yy y_right(:,2:end)];
    else
	yy = [yy(:,1:end-1) y_right];
    end
end

% Solve for the coefficients
x = yy\y;
if norm(yy*x-y) > 1e-4
    warning('Pulse response fit is poor.');
end
switch form 
    case 'FB'
	if ~iscell(tdac)
	    Bc2 = [ x(1:order) zeros(order,n_extra) ];
	    if (n_extra > 0)
		Dc2 = [ 0 x(order+1:end).'];
	    else
		Dc2 = x(order+1:end).';
	    end
	else
	    BcDc = [Bc;Dc];
	    i = find(BcDc);
	    BcDc(i) = x;
	    Bc2 = BcDc(1:end-1,:);
	    Dc2 = BcDc(end,:);
	end
    case 'FF'
	Bc2 = [Bc zeros(order,n_extra)];
	Cc = x(1:order).';
	if (n_extra > 0)
	    Dc2 = [ 0 x(order+1:end).'];
	else
	    Dc2 = x(order+1:end).';
	end
    otherwise
	fprintf(1,'%s error. No code for form "%s".\n', mfilename, form);
end
Dc1 = 0;
Dc = [Dc1 Dc2];
Bc1 = [1; zeros(order-1,1)];
Bc = [Bc1 Bc2];
% Scale Bc1 for unity STF magnitude at f0
fz = angle(ntf.z{1})/(2*pi);
f1 = fz(1);
ibz = abs(fz-f1) <= abs(fz+f1);
fz = fz(ibz);
f0 = mean(fz);
if min(abs(fz)) < 3*min(abs(fz-f0))
    f0 = 0;
end
L0c = zpk(ss(Ac,Bc1,Cc,Dc1));
G0 = evalTFP(L0c,ntf,f0);
if f0 == 0
    Bc(:,1) = Bc(:,1)*abs(Bc(1,2:end)*(tdac2(2:end,2)-tdac2(2:end,1))/Bc(1,1));
else
    Bc(:,1) = Bc(:,1)/abs(G0);
end

ABCDc = [Ac Bc; Cc Dc];
ABCDc = ABCDc .* ( abs(ABCDc) > eps^(1/2) ) ;

