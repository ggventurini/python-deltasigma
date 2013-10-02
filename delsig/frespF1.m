function fresp = frespF1(f1,f,phi)
%frespF1(f1,f)	Plot/calculate the frequency response of the F1 filter 
%in a Saramaki HBF at the points given in the optional f (n by 1) vector.
if nargin < 3
    phi = 1;
    if nargin < 2
	f = linspace(0,0.5);
    end
end

cos_w = cos(2*pi*f)';
F1 = 0.5;
for i = 1:length(f1);
    F1 = F1 + f1(i)*((cos_w/phi).^(2*i-1));
end

if nargout == 0
    plot(f,dbv(F1))
    grid on
else
    fresp = F1;
end
