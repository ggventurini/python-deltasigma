function un=uvar(u,N);
%un = uvar(u,N);	Compute a bounded sequence that has 90% of its 
%values at the extremes.

u1 = u(1); u2=u(2);
un = u1(ones(1,N));
r = rand(1,N);
ri2 = r > 0.55;
ri1 = r > 0.45 & ~ri2;
un(ri1) = u1 + (u2-u1)*10*(r(ri1)-.45);
un(ri2) = u2(ones(1,sum(ri2)));
