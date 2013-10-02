function s = find2dPIS(u,ABCD,options)
%function s = find2dPIS(u,ABCD,options)
%Find a positively invariant set for the 2nd-order binary modulator whose 
%loop filter is described by ABCD and whose input is a constant u.
%Return s = [] if no invariant set is found.
%
%options = [ dbg=0 itnLimit=100 expFactor=0.005 N=1000 skip=100]
%dbg=1 causes debugging information to be displayed.

s=[];
n = size(ABCD,1)-1;
if n ~= 2
    fprintf('%s: Error. The modulator must be second order.\n', mfilename);
    return;
end
% Handle the options
nlev=2;
parameters = ['dbg      ';'itnLimit ';'expFactor';'N        ';'skip     '];
defaults = [ 0 100 .005 1000 100];
for i=1:length(defaults)
    if i>length(options)
    	eval([parameters(i,:) '=defaults(i);'])  
    elseif isnan(options(i)) 
    	eval([parameters(i,:) '=defaults(i);'])  
    else
    	eval([parameters(i,:) '=options(i);'])  
    end
end

% Compute a few iterations of difference equations
if size(u)==[1 1]
    un = u(ones(1,N+skip));
elseif size(u) == [2 1]	
    if ABCD(n+1,n+1) ~= 0	% Require D1=0
    	fprintf('%s: Limitation. D1 must be zero for u-ranges.\n', mfilename);
    	return;
    end
    % Make 90% of the u-values take on the extremes
    un = uvar(u,skip+N);
else
    fprintf('%s: Error. Argument 2 (u) has the wrong dimensions.\n', mfilename);
    return
end
[v x] = simulateDSM(un,ABCD,nlev);
x = x(:,skip+1:skip+N);

xmin = min(x(1,:)); xmax = max(x(1,:)); dx = xmax-xmin;
ymin = min(x(2,:)); ymax = max(x(2,:)); dy = ymax-ymin;
axis1 = [ xmin-dx/4 xmax+dx/4 ymin-dy/4 ymax+dy/4 ];

% Take the convex hull of the result
s = hull2d(x')';
ec = mean(x')';

for i = 1:itnLimit
    % Inflate the hull
    shift = ec(:,ones(1,length(s)));
    s = shift + (1+expFactor)*( s-shift );

    % Split the set
    [splus eplus sminus eminus] = dssplit2d(u,ABCD,s);
    % Map the two halves 
    s1 = dsmap(u,ABCD,2,splus,eplus,1);
    s2 = dsmap(u,ABCD,2,sminus,eminus,-1);
    ns = [s1(:,1:size(s1,2)-1) s2(:,1:size(s2,2)-1)];

    % Test for inclusion: ns inside s (the inflated hull)
    out = outconvex2d(ns,s);
    if dbg
    	clf; hold on; grid; axis(axis1);
        set(gca,'drawmode','fast');
    	dotplot(x,'k.');
    	dotplot(ec,'o');
    	polyplot(s);
    	polyplot(s1,'m');
    	polyplot(s2,'c');
    	outi = logical(sign(out));
    	dotplot(ns(:,outi),'rs');
    	str = sprintf('Iteration %d: %d image vertices outside',i, sum(outi));
    	title(str);
        drawnow;
    end
    if out == 0
    	break;
    end
    % take the hull and repeat
    s = hull2d(ns')';
end
