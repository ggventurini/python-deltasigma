== Versions ==

The feature support, split by release, is expected as follows:

**V. 0.1** Support for all base band Delta-Sigma topologies. Continuous time and discrete time.

**V. 0.2** Support for Quadrature modulators.

**V. 0.3** PIS set calculation.

**V. 0.4** HBF and LC topologies.  

== Porting ==

The conversion will be performed according to the following steps:

1. **MATLAB/Python conversion**
  1. Convert the M-files to Python, whenever possible, to have several 
Python modules providing a 1:1 replacement of the original MATLAB
functions. 
  1. Whenever possible keep the exact function signatures, 
return values and default parameters. 
  1. Have a working unit test developed at the same time as each function.
  1. Have documentation developed at the same time as each function.
  1. Do not address the C Matlab EXtension (C MEX) files at this stage.

2. **C MEX / Cython conversion and other Cython files**
  1. Study the C MEX files. Can their functionality be reproduced in
straight Python, Python + numpy or Python + numpy/scipy + 
python-controls, without a *very* high penality in terms of 
speed and/or complexity? If yes, implement the Python module.
  2. For those modules which really need to be implemented in C as 
a Python extension, look into a possible Cython implemetation.
  3. If Cython is not an option, then consider Pyfort. 

3. **Drop dependencies, if possible**
  1. Do we really need all the current dependencies? Could we 
live without python-control? If yes, rewrite and remove the
superfluous modules.

5. Have fun.
