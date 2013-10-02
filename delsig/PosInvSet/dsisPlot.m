function dsisPlot(dbg,itn,order,x,s,e,ns,out)
% Show some pretty pictures
if dbg==0
    return;
elseif dbg==2
    fprintf('Iteration %d: %d/%d image vertices outside\n',itn, sum(out),size(ns,2));
    return
end

if itn>0
    str = sprintf('Iteration %d: %d/%d image points outside',itn, sum(out), size(ns,2));
else
    str = sprintf('Final object: %d vertices and %d edges', size(s,2), size(e,2));
end
clf;

if order == 2
    dotplot(x);
    hold on; grid;
    polyplot(s,'k');
    dotplot(ns,'ko');
    dotplot(ns(:,logical(out)),'rs');
    title(str);
    drawnow;
elseif order >= 3 
    dotplot(x);
    edgeplot(e,s);
    dotplot(ns,'+');
    dotplot(ns(:,logical(out)),'rs');
    subplot(222);
    title(str);
    drawnow;
end
