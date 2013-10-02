function ntf = synthesizeNTF0(order,OSR,opt,H_inf,f0)


% Determine the zeros.
if f0~=0		% Bandpass design-- halve the order temporarily.
    order = order/2;
    dw = pi/(2*OSR);
else
    dw = pi/OSR;
end

if length(opt)==1
    if opt==0
	z = zeros(order,1);
    else
	z = dw*ds_optzeros(order,opt);
	if isempty(z)
		    return;
	end
    end
    if f0~=0		% Bandpass design-- shift and replicate the zeros.
	order = order*2;
	z = z + 2*pi*f0;
	ztmp = [ z'; -z' ];
	z = ztmp(:);
    end
    z = exp(j*z);
else
    z = opt(:);
end
		
ntf = zpk(z,zeros(1,order),1,1);
itn_limit = 100;

% Iteratively determine the poles by finding the value of the x-parameter
% which results in the desired H_inf.
if f0 == 0			% Lowpass design
    HinfLimit = 2^order;  % !!! The limit is actually lower for opt=1 and low OSR
    if H_inf >= HinfLimit 
	fprintf(2,'%s warning: Unable to achieve specified Hinf.\n', mfilename);
	fprintf(2,'Setting all NTF poles to zero.\n');
	ntf.p = zeros(order,1);
    else
	x=0.3^(order-1);	% starting guess
	converged = 0;
	for itn=1:itn_limit
	    me2 = -0.5*(x^(2./order));
	    w = (2*[1:order]'-1)*pi/order;
	    mb2 = 1+me2*exp(j*w);
	    p = mb2 - sqrt(mb2.^2-1);
	    out = find(abs(p)>1);
	    p(out) = 1./p(out);	% reflect poles to be inside the unit circle.
	    ntf.z = z;	ntf.p = cplxpair(p);
	    f = real(evalTF(ntf,-1))-H_inf;
	% [ x f ]
	    if itn==1 
		delta_x = -f/100;
	    else
		delta_x = -f*delta_x/(f-fprev);
	    end
	    
	    xplus = x+delta_x;
	    if xplus>0 
		x = xplus;
	    else
		x = x*0.1;
	    end
	    fprev = f;

	    if abs(f)<1e-10 | abs(delta_x)<1e-10 
		converged = 1;
		break;
	    end
	    if x>1e6
		fprintf(2,'%s warning: Unable to achieve specified Hinf.\n', mfilename);
		fprintf(2,'Setting all NTF poles to zero.\n');
		ntf.z = z;	ntf.p = zeros(order,1);
		break;
	    end
	    if itn == itn_limit
		fprintf(2,'%s warning: Danger! Iteration limit exceeded.\n',...
	 	  mfilename);
	    end
	end
    end
else				% Bandpass design.
    x = 0.3^(order/2-1);	% starting guess (not very good for f0~0)
    if f0>0.25
	z_inf=1;
    else
	z_inf=-1;
    end
    c2pif0 = cos(2*pi*f0);
    for itn=1:itn_limit
	e2 = 0.5*x^(2./order);
	w = (2*[1:order]'-1)*pi/order;
	mb2 = c2pif0 + e2*exp(j*w);
	p = mb2 - sqrt(mb2.^2-1);
	% reflect poles to be inside the unit circle.
	out = find(abs(p)>1);
	p(out) = 1./p(out);
	ntf.z = z;	ntf.p = cplxpair(p);
	f = real(evalTF(ntf,z_inf))-H_inf;
% 	[x f]
	if itn==1 
	    delta_x = -f/100;
	else
	    delta_x = -f*delta_x/(f-fprev);
	end
	
	xplus = x+delta_x;
	if xplus > 0
	    x = xplus;
	else
	    x = x*0.1;
	end

	fprev = f;
	if abs(f)<1e-10 | abs(delta_x)<1e-10
	    break;
	end
	if x>1e6
	    fprintf(2,'%s warning: Unable to achieve specified Hinf.\n', mfilename);
	    fprintf(2,'Setting all NTF poles to zero.\n');
	    ntf.p = zeros(order,1);
	    break;
	end
	if itn == itn_limit
	    fprintf(2,'%s warning: Danger! Iteration limit exceeded.\n',...
	      mfilename);
	end
    end
end

z = cplxpair(ntf.z{:});
p = cplxpair(ntf.p{:});
rev = order:-1:1;

% Assemble the ntf struct
ntf.z = z(rev)';
ntf.p = p(rev)';
