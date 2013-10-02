function edgeplot(e,s,fmt)
%function edgeplot(e,s,fmt) Plot edges e with format fmt
if nargin < 3 
    fmt = '-';
end
if size(s,1) == 2
    hold on; grid on; 
    % Assume that the s points are in order.
    plot([s(1,:) s(1,1)], [s(2,:) s(2,1)], fmt);
elseif size(s,1) == 3
    tbl = [ 1 1 2; 4 1 3; 3 2 3; 2 9 9 ];
    for p=1:3
	subplot(2,2,tbl(p,1))
	hold on; grid on; 
	x = tbl(p,2); y = tbl(p,3);
	xlabel(['x' num2str(x)]); ylabel(['x' num2str(y)]);
	for i=1:size(e,2)
	    p1 = s(:,e(1,i)); p2 = s(:,e(2,i));
	    plot([p1(x) p2(x)], [p1(y) p2(y)], fmt);
	end
    end
    subplot(2,2,tbl(4,1));
    hold on; grid on;
    xlabel('x1'); ylabel('x2'); zlabel('x3');
    for i=1:size(e,2)
	p1 = s(:,e(1,i)); p2 = s(:,e(2,i));
	plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], fmt);
    end
elseif size(s,1) == 4
    tbl = [ 1 1 2 3; 2 1 2 4; 3 1 3 4; 4 2 3 4 ];
    for p=1:4
	subplot(2,2,tbl(p,1))
	hold on; grid on; 
	x = tbl(p,2); y = tbl(p,3); z = tbl(p,4);
	xlabel(['x' num2str(x)]); ylabel(['x' num2str(y)]); zlabel(['x' num2str(z)]);
	for i=1:size(e,2)
	    p1 = s(:,e(1,i)); p2 = s(:,e(2,i)); 
	    plot3([p1(x) p2(x)], [p1(y) p2(y)], [p1(z) p2(z)], fmt);
	end
    end
end
