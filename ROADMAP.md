## Versions

The feature support, split by release, is expected as follows:

**V. 0.1** Support for all base band Delta Sigma topologies. Continuous time and discrete time.

**V. 0.2** Support for Quadrature modulators.

**V. 0.4** LC topologies.

**V. 0.3** HBF and PIS set calculation.

## Porting

The conversion will be performed according to the following steps:

1. **MATLAB/Python conversion**
  1. Convert the M-files to Python, whenever possible, to have several Python
     modules providing a 1:1 replacement of the original MATLAB functions.
  1. Whenever possible keep the exact function signatures, return values and
     default parameters.
  1. Have a working unit test developed at the same time as each function.
  1. Have documentation developed at the same time as each function.
  1. Do not address the C Matlab EXtension (C MEX) files at this stage.

2. **C MEX / Cython conversion and other Cython files**
  2. For those modules which need to be implemented in C as a Cython extension,
     look into a possible Cython implementation.

3. Have fun.
