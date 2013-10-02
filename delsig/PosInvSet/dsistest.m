% parameters which may be defaulted: u, order, dbg, expFactor
if exist('u')~=1,		u = 1/pi;		end
if exist('order')~=1,		order = 2;		end
if exist('dbg')~=1,		dbg = 1;		end
if exist('expFactor')~=1	expFactor=0.01;		end
if exist('nlev')~=1,		nlev = 2;		end
if exist('qhullArgA')~=1,	qhullArgA=0.999;	end
if exist('qhullArgC')~=1,	qhullArgC=.001;		end

figure(1);
if order==2
    mod = mod2;
    ABCD = mod.ABCD;
    t = cputime;
    [s e n o] = findPIS(u,ABCD,nlev, ...
	[dbg NaN expFactor NaN NaN qhullArgA qhullArgC]);
    t = cputime - t
else 
    H = synthesizeNTF(order);
    [a,g,b,c] = realizeNTF(H);
    ABCD = stuffABCD(a,g,b,c);
    t = cputime;
    [s e n o] = findPIS(u,ABCD,nlev, ...
	[dbg NaN expFactor NaN NaN qhullArgA qhullArgC]);
    t = cputime - t
end

figure(2);
clf;
if isempty(n)
    fprintf('No positively invariant set found.\n');
else
    N=10000; skip=100;
    if length(u)==1
    	uu = u(ones(1,N+skip));
    else
    	uu = uvar(u,N+skip);
    end
    [junk x] = simulateDSM(uu,ABCD,nlev);
    x = x(:,1+skip:N+skip);

    ns = dsmap(u,ABCD,nlev,s,e);        
    out = outsideConvex(ns,n,o);         
    dsisPlot(1,0,order,x,s,e,ns,out);
    if(order==3)
    	subplot(222);
    	view(20,20);
    end;

    fprintf('%d points from the %d simulated states are outside.\n', sum(outsideConvex(x,n,o)),N);         
    fprintf('%d image points are outside.\n', sum(out));         
    fprintf('The returned object has %d vertices, %d edges and %d faces.\n', ...
      size(s,2),size(e,2),size(n,2));
    %fprintf('order=%d, u=%f, expFactor=%f\n', order, u, expFactor );
end
