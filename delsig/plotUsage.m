function plotUsave(sv)
% Plot the elemet usage for a multi-elemet DAC
T = size(sv,2);					M = size(sv,1);

% Plot the grid
x=[0:T; 0:T];					x=x(:)';
T2 = ceil((T+1)/2);
y=[zeros(1,T2); M*ones(2,T2); zeros(1,T2)];	y = y(1:2*(T+1))';
plot(x,y,'k');
hold on;
M2 = ceil((M+1)/2);
x=[zeros(1,M2); T*ones(2,M2); zeros(1,M2)];	x = x(1:2*(M+1))';
y=[0:M; 0:M];					y=y(:)';
plot(x,y,'k');
% axis([0 T 0 M]);
axis('image');

for t=1:T
    for i=1:M
	if sv(i,t) == 1
	    fill([t-1 t-1 t t],[i-1 i i i-1],'b');
	elseif sv(i,t) == -1
	    fill([t-1 t-1 t t],[i-1 i i i-1],'g');
	elseif sv(i,t) == 1i
	    fill([t-1 t-1 t t],[i-1 i i i-1],'r');
	elseif sv(i,t) == -1i
	    fill([t-1 t-1 t t],[i-1 i i i-1],'y');
	end
    end
end
hold off
% xlabel('time');
% ylabel('element number');
