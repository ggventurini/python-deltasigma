% demonstrate findPIS
addPIS;		% Add the PosInvSet directory to MATLAB's path
clear;
if exist('qhull')==3
    clc;
    fprintf(1,'\t\tInvariant Set for MOD2 (findPIS) \n\n');
    order=2;
    dsistest
    fprintf(1,'paused\n');
    pause

    clc;
    fprintf(1,'\t\tInvariant Set for a 3rd-order modulator (findPIS) \n\n');
    order=3;
    dsistest
else
    clc;
    fprintf(1,'\t\tInvariant Set for MOD2 (find2dPIS) \n\n');
    mod = mod2;
    ABCD = mod.ABCD; [A B C D] = partitionABCD(ABCD);
    order=2;
    figure(1);
    set(gcf,'NumberTitle','off'); 
    set(gcf,'Name','Working Invariant Set');
    echo on;
    u = 1/pi;
    t = cputime;
    s = find2dPIS(u,ABCD,1)
    t = cputime-t
    echo off;
    N=10000; skip=100;
    [junk x] = simulateDSM(u(ones(1,N+skip)),ABCD,2);
    x = x(:,1+skip:N+skip);
    nv = size(s,2);
    [splus, eplus, sminus, eminus] = dssplit2d(u,ABCD,s);
    Buv = B*[u;1];
    s1 = A*splus + Buv(:,ones(1,size(splus,2)));        
    Buv = B*[u;-1];
    s2 = A*sminus + Buv(:,ones(1,size(sminus,2)));       
    ns = [s1 s2];                        
    out = outconvex2d(ns,s);

    figure(2);
    set(gcf,'NumberTitle','off'); 
    set(gcf,'Name','Final Invariant Set');
    clf; hold on; grid;
    dotplot(x,'k.');
    polyplot(s,'b');
    polyplot(s1,'m');
    polyplot(s2,'c');
    outi = logical(sign(out));
    dotplot(ns(:,outi),'rs');
    str = sprintf('Final Object: %d image vertices outside', sum(outi));
    title(str);
    xlabel('x_1');
    ylabel('x_2');
    drawnow;

    fprintf(1,'%d points from the %d simulated states are outside.\n', sum(outconvex2d(x,s)),N);
    fprintf(1,'%d image points are outside.\n', sum(out));         
    fprintf(1,'The returned polygon has %d vertices.\n', size(s,2));
end

