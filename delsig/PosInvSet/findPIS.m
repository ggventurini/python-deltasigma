function [s,e,n,o,Sc] = findPIS(u,ABCD,nlev,options)
%[s e n o Sc] = findPIS(u,ABCD,nlev=2,options). Find a positively-invariant set.
%findPIS finds a convex positively invariant set for the (nlev)-level
%delta-sigma modulator described by ABCD whose input is u.
%u may be either a scalar or a 2x1 vector. 
%In the former case the input is constant;
%in the latter the input is a sequence of arbitrary values in the given range.
%The invariant set is described with the following parameters:
%s = its vertices, e = its edges, n = facet normals, o = facet offsets. 
%Points inside the set are characterized by  n'*x + o <= 0.
%
% options = [ dbg=0 itnLimit=2000 expFactor=0.01 N=1000 skip=100
%	      qhullArgA=0.999 qhullArgC=.001]
%qhullArgA is the cosine of the angle tolerance (controls facet-merging).
%qhullArgC is the centrum distance tolerance (controls facet-merging).
%dbg=1 causes debugging information to be displayed.
if nargin>2
    if isnan(nlev) | isempty(nlev)
	nlev = 2;
    end
end
order = size(ABCD,1)-1;

% Handle the options
optionName = ['dbg      ';'itnLimit ';'expFactor';'N        ';'skip     ';'qhullArgA';'qhullArgC'];
defaults = [ 0 2000 .01 1000 100 0.999 0.001];
for i=1:length(defaults)
    if i>length(options)
	eval([optionName(i,:) '=defaults(i);'])  
    elseif isnan(options(i)) 
	eval([optionName(i,:) '=defaults(i);'])  
    else
	eval([optionName(i,:) '=options(i);'])  
    end
end
qhullArgs = sprintf('qhull Qcx A%g C%g', qhullArgA, qhullArgC); 

% Compute a few iterations of difference equations
if size(u)==[1 1]
    un = u(ones(1,skip+N));
elseif size(u) == [2 1]	
    if ABCD(order+1,order+1) ~= 0	% Require D1=0
	fprintf('%s: Limitation. D1 must be zero for u-ranges.\n', mfilename);
	return;
    end
    % Make 90% of the u-values take on the extremes
    un = uvar(u,skip+N);
else
    fprintf('%s: Error. Argument 1 (u) has the wrong dimensions.\n', mfilename);
    return
end
[v x xmax] = simulateDSM(un,ABCD,nlev);
if max(xmax) > 100
    fprintf('%s: A direct simulation indicates that the modulator is unstable.\n', mfilename);
    s = Inf;
    s = s(ones(order,1));
    e=[]; n=[]; o=[];
    return
end
x = x(:,1+skip:N+skip);

%Do scaling (coordinate transformation) to help qhull do better facet merging.
%The scaling is based on principal component analysis (pg 105 of notebook 6).
center = mean(x')';
xp = x-center(:,ones(1,size(x,2)));
R = xp*xp'/N;
[Q L] = eig(R);
Sc = Q*sqrt(L); Si = inv(Sc);
[A0 B0 C0 D0] = partitionABCD(ABCD);
ABCD = [ Si*A0*Sc Si*B0; C0*Sc D0];
x = Si*x;
xmax = max(abs(x)')';
center = Si*center;

%Store original data in case I need to restart
restart = 1;
x0 = x;
ABCD0 = ABCD;
Si0 = Si; Sc0 = Sc;
xmax0 = xmax; center0 = center;

converged = 0;
for itn = 1:itnLimit
    if restart==1
	restart = 0;
	% Use the hull of x for the first iteration.
	[s e n o] = qhull(x,qhullArgs);
    else
	% Expand the outside points
	ns = dsexpand(ns(:,logical(out)),center,expFactor);
	% Use the hull of s and the expanded ns for the next iteration.
	[s e n o] = qhull([s ns], qhullArgs);
    end
    % Map the set
    ns = dsmap(u,ABCD,nlev,s,e);
    if ~isempty(find(max(abs(ns'))' > 10*xmax))
	fprintf('Set is much larger than necessary--\n');
	fprintf('Reducing expansion factor, increasing hull accuracy, and re-starting.\n');
	restart = 1;
	expFactor = 0.5*expFactor; 
	qhullArgC = 0.5*qhullArgC;
	qhullArgA = 0.75 + 0.25*qhullArgA;
	qhullArgs = sprintf('qhull Qcx A%g C%g', qhullArgA, qhullArgC); 
	x = x0; xmax = xmax0; center = center0;
	Si = Si0; Sc = Sc0; ABCD = ABCD0;
    else
	% Test for inclusion: ns inside s 
	out = outsideConvex(ns,n,o);
	% Draw some pretty pictures or print some status information.
	dsisPlot(dbg,itn,order,x,s,e,ns,out);
	if out == 0
	    % Check the PIS by forming the exact hull 
	    % and checking its images for inclusion.
	    [ss ee nn oo] = qhull(s,'exact Qcx');
	    ns = dsmap(u,ABCD,nlev,ss,ee);
	    out =  outsideConvex(ns,nn,oo);
	    if( sum(out) == 0 )
		converged = 1;
		break;
	    end
	    fprintf('Apparent convergence, but %d vertices outside.\n',sum(out));
	    fprintf('Continuing with tighter hull tolerances.\n');
	    % Halve the centrum distance and inter-normal angle parameters.
	    qhullArgC = 0.5*qhullArgC;
	    qhullArgA = 0.75 + 0.25*qhullArgA;
	    qhullArgs = sprintf('qhull Qcx A%g C%g', qhullArgA, qhullArgC); 
	end

	center = mean(ns')';
	% The following can be done intermittently
	N = size(ns,2);
	xp = ns-center(:,ones(1,N));
	R = xp*xp'/N; [Q L] = eig(R);
	% Re-do the scaling if the principal axis is too long.
	if max(max(L)) > 1.5	
	    if dbg>1
		fprintf('Re-doing the scaling at iteration %d.\n', itn);
	    end
	    Sc1 = Q*sqrt(L); Si1 = inv(Sc1);
	    Sc = Sc*Sc1; Si = inv(Sc);
	    s = Si1*s; ns = Si1*ns; x = Si1*x; center = Si1*center;
	    xmax = max(abs(x)')';	% !! I should use the hull of the points
	    ABCD = [ Si*A0*Sc Si*B0; C0*Sc D0];
	end
    end
end	% for itn...

if converged
    % Undo the scaling.
    s = Sc*s;
    n = Si'*n;
    Sn = 1 ./ sqrt(sum(n.^2)); 
    % n = n * diag(Sn); This is inefficient if the number of faces is large.
    for i=1:size(n,2)
        n(:,i) = n(:,i)*Sn(i);
    end
    o = o .* Sn;
else
    fprintf('%s: Unable to determine stability.\n', mfilename);
    s = Inf;
    s = s(ones(order,1));
    e=[]; n=[]; o=[];
end
