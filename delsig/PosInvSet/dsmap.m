function nx = dsmap(u,ABCD,nlev,x,e,v)
%function nx = dsmap(u,ABCD,nlev,x,e,v)
%For a DSM with input u, a structure ABCD and an nlev-level quantizer,
%compute the (potential) vertices of the image of a convex object 
%described in terms of its vertices (x) and edges (e).
%If u has two elements, it is considered to represent a range.
%v is the assumed quantizer output; it is computed if it is not supplied.

%Basic Algorithm:
% 1) Compute the images of the vertices.
% 2) For those edges which cross a splitting hyperplane,
%    compute the images of each split point.
% 3) For u-ranges, append the translated images to the list.

n = size(ABCD,1)-1;
A = ABCD(1:n, 1:n);
B = ABCD(1:n, n+1:n+2);
C = ABCD(n+1, 1:n);
D1= ABCD(n+1, n+1);	% For any DSM, D2=ABCD(n+1,n+2) must be zero.

N = size(x,2);
if length(u)==2
    u2 = u(2);
    u  = u(1);
    isRange = 1;
    if D1 ~=0
    	fprintf('%s: Limitation: D1 must be zero.\n');
	return;
    end
elseif length(u)==1
    isRange = 0;
else
    fprintf('%s: Error. The dimensions of u are wrong.\n');
    return;
end

% Compute v. The assumption that D1=0 for u-ranges is implicit in this step.
if nargin < 6
    y = C*x + D1*u;
    v = ds_quantize(y,nlev);
elseif length(v)~=N
    v = v(ones(1,N));
else
    fprintf('%s error: the supplied v argument is the wrong size.\n',mfilename);
end

% 1) Compute the images of the vertices.
B1u = B(:,1)*u;
nx = A*x + B1u(:,ones(1,N)) + B(:,2)*v;

% 2) For those edges which cross a (or several) splitting hyperplanes,
%    compute the two images of each split point.
diff = abs(v(e(1,:))-v(e(2,:)));
split1 = diff==2;	%edges split in one place only

% Handle the split1 edges en masse.
if any(split1)
    i1 = e(1,split1); i2 = e(2,split1);
    y0 = 0.5*(v(i1)+v(i2));	%the appropriate quantizer thresholds
    k1 = (y(i2)-y0)./(y(i2)-y(i1)); k2 = 1-k1;
    psplit = k1(ones(n,1),:).*x(:,i1) + k2(ones(n,1),:).*x(:,i2);
    N = length(k1);
    images1 = A*psplit + B1u(:,ones(1,N)) + B(:,2)*v(i1);
    images2 = images1 + B(:,2)*(v(i2)-v(i1));
    nx = [nx images1 images2];
end

% Treat the multiply-split edges as a special case.
split2 = find(diff>2);
for i = split2
    i1 = e(1,i);	i2 = e(2,i);
    v1 = v(i1);		v2 = v(i2);
    x1 = x(:,i1);	x2 = x(:,i2);
    y1 = y(i1);		y2 = y(i2);
    dv = v2-v1;		N = abs(dv/2);
    y0 = v1+sgn(dv);			%the first quantizer threshold crossed
    k1 = (y2-y0)/(y2-y1);	k2 = 1-k1;
    x0 = k1*x1 + k2*x2;			%the first split point
    image0 = A*x0 + B1u + B(:,2)*v1;	%its image
    deltaB = B(:,2)*(2*dv);		%the image shift due to splitting
    A_deltax = A*((x2-x1)/(0.5*(y2-y1)));	% The image shift due to x
    images = image0(:,ones(1,N)) + A_deltax*[0:N-1]; 
    images = [images images+deltaB(:,ones(1,N))];
end

% 3) For u-ranges, append the translated images to the list.
if isRange
    translation = (u2-u)*ABCD(1:n,n+1);
    nx = [nx nx+translation( :, ones(1,size(nx,2)) )];
end

