function ntf = synthesizeNTF1(order,osr,opt,H_inf,f0)


% Determine the zeros.
if f0~=0		% Bandpass design-- halve the order temporarily.
    order = order/2;
    dw = pi/(2*osr);
else
    dw = pi/osr;
end

if length(opt)==1
    if opt==0
        z = zeros(order,1);
    else
        z = dw*ds_optzeros(order,1+rem(opt-1,2));
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

zp = z(angle(z)>0);
x0 = (angle(zp)-2*pi*f0) * osr / pi;
if opt==4 & f0~=0
    % Do not optimize the zeros at f0
    x0(find( abs(x0)<1e-10 )) = [];
end

ntf = zpk(z,zeros(1,order),1,1);
Hinf_itn_limit = 100;

opt_iteration = 5;	% Max number of zero-optimizing/Hinf iterations
while opt_iteration > 0
    % Iteratively determine the poles by finding the value of the x-parameter
    % which results in the desired H_inf.
    ftol = 1e-10;
    if f0>0.25
        z_inf=1;
    else
        z_inf=-1;
    end
    if f0 == 0			% Lowpass design
        HinfLimit = 2^order;  % !!! The limit is actually lower for opt=1 and low osr
        if H_inf >= HinfLimit
            fprintf(2,'%s warning: Unable to achieve specified Hinf.\n', mfilename);
            fprintf(2,'Setting all NTF poles to zero.\n');
            ntf.p = zeros(order,1);
        else
            x=0.3^(order-1);	% starting guess
            converged = 0;
            for itn=1:Hinf_itn_limit
                me2 = -0.5*(x^(2./order));
                w = (2*[1:order]'-1)*pi/order;
                mb2 = 1+me2*exp(j*w);
                p = mb2 - sqrt(mb2.^2-1);
                out = find(abs(p)>1);
                p(out) = 1./p(out);	% reflect poles to be inside the unit circle.
                p = cplxpair(p);
                ntf.z = z;	ntf.p = p;
                f = real(evalTF(ntf,z_inf))-H_inf;
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
                if abs(f)<ftol | abs(delta_x)<1e-10
                    converged = 1;
                    break;
                end
                if x>1e6
                    fprintf(2,'%s warning: Unable to achieve specified Hinf.\n', mfilename);
                    fprintf(2,'Setting all NTF poles to zero.\n');
                    ntf.z = z;	ntf.p = zeros(order,1);
                    break;
                end
                if itn == Hinf_itn_limit
                    fprintf(2,'%s warning: Danger! Iteration limit exceeded.\n',...
                        mfilename);
                end
            end
        end
    else				% Bandpass design.
        x = 0.3^(order/2-1);	% starting guess (not very good for f0~0)
        c2pif0 = cos(2*pi*f0);
        for itn=1:Hinf_itn_limit
            e2 = 0.5*x^(2./order);
            w = (2*[1:order]'-1)*pi/order;
            mb2 = c2pif0 + e2*exp(j*w);
            p = mb2 - sqrt(mb2.^2-1);
            % reflect poles to be inside the unit circle.
            out = find(abs(p)>1);
            p(out) = 1./p(out);
            p = cplxpair(p);
            ntf.z = z;	ntf.p = p;
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
            if abs(f)<ftol | abs(delta_x)<1e-10
                break;
            end
            if x>1e6
                fprintf(2,'%s warning: Unable to achieve specified Hinf.\n', mfilename);
                fprintf(2,'Setting all NTF poles to zero.\n');
                p = zeros(order,1);
                ntf.p = p;
                break;
            end
            if itn == Hinf_itn_limit
                fprintf(2,'%s warning: Danger! Hinf iteration limit exceeded.\n',...
                    mfilename);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opt < 3   % Do not optimize the zeros
        opt_iteration = 0;
    else
        if f0 == 0
            ub = ones(size(x0));
            lb = zeros(size(x0));
        else
            ub = 0.5*ones(size(x0));
            lb = -ub;
        end
        options = optimset('TolX',0.001, 'TolFun',0.01, 'MaxIter',100 );
        options = optimset(options,'LargeScale','off');
        options = optimset(options,'Algorithm','active-set');
        options = optimset(options,'Display','off');
        %options = optimset(options,'Display','iter');
        vn = sscanf(version,'%d');
        if vn>=6
            x = fmincon( @(x) ds_synNTFobj1(x,p,osr,f0), x0,[],[],[],[], lb,ub,[],options);
        else
            error('To use opt>=3 you need to have the Optimization Toolbox and Matlab 6 or higher');
        end
        x0 = x;
        z = exp(2i*pi*(f0+0.5/osr*x));
        if f0>0
            z = padt(z,length(p)/2,exp(2i*pi*f0));
        end
        z = [z conj(z)]; z = z(:);
        if f0==0
            z = padt(z,length(p),1);
        end
        ntf.z = z;	ntf.p = p;
        if  abs( real(evalTF(ntf,z_inf)) - H_inf ) < ftol
            opt_iteration = 0;
        else
            opt_iteration = opt_iteration - 1;
        end
    end
end
