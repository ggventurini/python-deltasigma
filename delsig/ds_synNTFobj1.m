function y = ds_synNTFobj1(x,p,osr,f0)
% y=ds_synNTFobj1(x,p,osr,f0)	Objective function for synthesizeNTF() 
z = exp(2i*pi*(f0+0.5/osr*x));
if f0>0
    z = padt(z,length(p)/2,exp(2i*pi*f0));
end 
z = [z conj(z)]; z = z(:);
if f0==0
    z = padb(z,length(p),1);
end
[f1,f2] = ds_f1f2(osr,f0);
ntf = zpk(z,p,1,1);
y = db(rmsGain(ntf,f1,f2));
