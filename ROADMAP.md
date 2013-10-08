The conversion will be performed according to the following steps:

1. Convert the M-files to Python, whenever possible to have working 
Python modules without touching the C Matlab EXtension (C MEX) files.
Have a working unit test for each file.

2. Arrange the converted files in a module (as opposed to the the 
current arrangement where each function has its own module.)

2. Study the C MEX files. Can their functionality be reproduced in
straight Python, Python + numpy or Python + numpy/scipy + 
python-controls, without a *very* high penality in terms of 
speed and/or complexity? If yes, implement the Python module.

3. For those modules which really need to be implemented in C as 
a Python extension, look into a possible Cython implemetation.

4. If Cython is not an option, then consider Pyfort. 

5. Have fun.
