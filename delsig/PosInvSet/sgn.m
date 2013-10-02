function y=sgn(x)
%y=sgn(x) The signum function. Unlike sign(), sgn(0)=1.
y=sign(sign(x)+.5);
