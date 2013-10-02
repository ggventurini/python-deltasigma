function F0 = evalF0(f1,z,phi)
%F0 = evalF0(f1,z,phi)	Calculate the values of the F0 (prototype) filter 
%of a Saramaki HBF at the given points.
F0 = evalF1( f1, 0.5*(z + 1./z), phi );
