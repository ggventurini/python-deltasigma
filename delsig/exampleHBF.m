function [F1, F2] = exampleHBF(n)
if nargin<1
    n = 1
end
switch n
    case 1
      [f1,q1,f2,q2] = HBFex1;
    case 2
      [f1,q1,f2,q2] = HBFex2;
    case 3
      [f1,q1,f2,q2] = Saramaki88;
    case 4
      [f1,q1,f2,q2] = Saramaki90;
    otherwise
      fprintf(2,'%s error: No code for example #%d available.\n',mfilename,n);
end
F1 = convertForm(f1,q1);
F2 = convertForm(f2,q2);
return

function F = convertForm(f,q)
n = length(f);
F = struct('val',cell(1,n),'csd',cell(1,n));
for i=1:n
    F(i).val = f(i);
    csdrep = q(2*i-1:2*i,:);
    F(i).csd = csdrep(:,csdrep(1,:)~=0);
end
return

function [f1,q1,f2,q2] = HBFex1
q1 = [	 0 -4 -5
	 1  1  1
	 0 -4 -5
	-1 -1 -1
	-1 -3 -5
	 1  1  1
	-3 -5  0
	-1 -1  0 ];
f1 = bunquantize(q1);
q2 = [	-1 -3 -7
	 1  1  1 
	-2 -5 -7
	-1  1  1
	-3 -8  0
	 1 -1  0
	-4 -6 -8
	-1 -1 -1
	-4 -8  0
	 1 -1  0
	-5 -7 -8
	-1 -1 -1
	-5  0  0
	 1  0  0
	-6 -7  0
	-1 -1  0
	-7 -8  0
	-1 -1  0
        -6  0  0
	 1  0  0 ];
f2 = bunquantize(q2);
return

function [f1,q1,f2,q2] = HBFex2
%Coefficients from toolbox v2.0 documentation
q1 = [	 0 -4 -7
	 1 -1  1
	-1 -3 -6
	-1 -1 -1
	-2 -4 -7
	 1 -1  1 ];
f1 = bunquantize(q1);
q2 = [	-1 -3 -8
	 1  1 -1
	-2 -4 -9
	-1  1 -1
	-3 -5 -9
	 1 -1  1
	-4 -7 -8
	-1  1  1
	-5 -8 -11
	 1 -1 -1
        -6 -9 -11
        -1  1 -1 ];
f2 = bunquantize(q2);
return

function [f1,q1,f2,q2] = Saramaki88()
% Coefficients from "Efficient VLSI-Realizable Decimators for Sigma-Delta
% Analog-to-Digital Converters," T. Saramaki and H. Tenhunen, ISCAS 1988,
% pp 1525-2528.
q1 = [	 0 -2
	 1 -1
	-2 -9
	-1  1 ];
f1 = bunquantize(q1);
q2 = [	-1 -3  0 
	 1  1  0 
	-3 -5 -6
	-1 -1 -1
	-4 -6  0
	 1  1  0 ];
f2 = bunquantize(q2);
return

function [f1,q1,f2,q2] = Saramaki90()
% Coefficients from "Multiplier-Free Decimation Algorithms for Superresolution
% Oversampled Converters," T. Saramaki, T. Karema, T. Ritoniemi, and H. Tenhunen
% ISCAS 1990, vol. 4. pp 1525-2528.
q1 = [	 0 -4 -5
	 1  1  1
	 0 -4 -5
	-1 -1 -1
	-1 -3 -5
	 1  1  1
	-3 -5 -20
	-1 -1  1 ];
f1 = bunquantize(q1);
q2 = [	-1 -3 -7
	 1  1  1 
	-2 -5 -7
	-1  1  1
	-3 -8  0
	 1 -1  0
	-4 -6 -8
	-1 -1 -1
	-4 -8  0
	 1 -1  0
	-5 -7 -8
	-1 -1 -1
	-5  0  0
	 1  0  0
	-6 -7  0
	-1 -1  0
	-6  0  0
	 1  0  0
	-7 -8  0
	-1 -1  0
	-6  0  0
	 1  0  0];
f2 = bunquantize(q2);
return
