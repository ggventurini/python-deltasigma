function optZeros = ds_optzeros( n, opt )
%optZeros = ds_optzeros( n, opt=1 )
%A helper function for the synthesizeNTF function of the Delta-Sigma Toolbox.
%Returns the zeros which minimize the in-band noise power of 
%a delta-sigma modulator's NTF.
if nargin<2
    opt = 1;
end

if opt==0
    optZeros = zeros(1,ceil(n/2));
else
    switch n
    case 1
	optZeros=0;
    case 2
	if opt==1
	    optZeros=sqrt(1/3);
	else
	    optZeros=0;
	end
    case 3
	optZeros=[sqrt(3/5) 0];
    case 4
	if opt==1
	    discr = sqrt(9./49-3./35);
	    tmp = 3./7;
	    optZeros=sqrt([tmp+discr tmp-discr]);
	else
	    optZeros=[0 sqrt(5/7)];
	end
    case 5
	discr = sqrt(25./81-5./21);
	tmp = 5./9;
	optZeros=sqrt([tmp+discr tmp-discr 0]);
    case 6
	if opt==1
	    optZeros=[ 0.23862059 0.66120988 0.9324696 ];
	else
	    discr = sqrt(56.)/33;
	    tmp = 7./11;
	    optZeros=sqrt([0 tmp+discr tmp-discr]);
	end
    case 7
	optZeros=[0 0.40584371 0.74153078 0.94910785 ];
    case 8
	if opt==1
	    optZeros=[ 0.18343709 0.52553345 0.79666684 0.96028993 ];
	else
	    optZeros=[ 0 0.50563161 0.79017286 0.95914731 ];
	end
    case 9
	optZeros=[ 0 0.32425101 0.61337056 0.83603082 0.9681602 ];
    case 10
	if opt==1
	    optZeros=[ 0.1834370913 0.5255334458 0.7966668433 0.9602899327];
	else
	    optZeros=[0 0.41572267 0.67208682 0.86238894 0.97342121 ];
	end
    case 11
	optZeros=[0 0.26953955 0.51909468 0.73015137 0.88706238 0.97822864];
    case 12
	if opt==1
	    optZeros=[0.12523875 0.36783403 0.58731921 0.7699033 0.90411753 0.9815607];
	else
	    optZeros=[0 0.35222363 0.58006251 0.76647993 0.90281326 0.98132047];
	end
    case 13
	optZeros=[ 0 0.23045331 0.44849063 0.64234828 0.8015776 0.91759824 0.98418306 ];
    case 14
	if opt==1
	    optZeros=[ 0.10806212 0.31911586 0.51525046 0.68729392 0.82720185 0.92843513 0.98628389 ];
	else
	    optZeros=[ 0 0.30524384 0.50836649 0.6836066 0.82537239 0.92772336 0.98615167 ];
	end
    otherwise
	fprintf(1,'Optimized zeros for n>14 are not available.\n');
	return;
    end
end

% Sort the zeros and replicate them.
z = sort(optZeros);
optZeros = zeros(n,1);
m=1;
if(rem(n,2)==1)
    optZeros(1) = z(1);
    z = z(2:length(z));
    m=m+1;
end
for(i=1:length(z))
    optZeros(m)   =  z(i);
    optZeros(m+1) = -z(i);
    m = m+2;
end

