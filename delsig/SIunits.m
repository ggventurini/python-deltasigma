function [factor,prefix] = SIunits(x)
% [factor,prefix] = SIunits(x)
% The factors and suffixes supported are:
% 1e-3	m	milli		1e3	k	kilo
% 1e-6	u	micro		1e6	M	mega
% 1e-9	n	nano		1e9	G	giga
% 1e-12	p	pico		1e12	T	tera
% 1e-15	f	femto		1e15	P	peta
% 1e-18	a	atto		1e18	E	exa
% 1e-21	z	zepto		1e21	Z	zeta
% 1e-24 y	yocto		1e24	Y	yotta
prefixes_n = {'m';'u';'n';'p';'f';'a';'z';'y'};
prefixes_p = {'k';'M';'G';'T';'P';'E';'Z';'Y'};
if x~=0
    p = floor( log10( abs(x) ) / 3 + eps);
else
    p = 0;
end
if p>0
    p = min(p,length(prefixes_p));
    prefix = prefixes_p{p};
    factor = 10^(3*p);
elseif p<0
    p = min(-p,length(prefixes_n));
    prefix = prefixes_n{p};
    factor = 10^(-3*p);
else
    prefix = '';
    factor = 1;
end

