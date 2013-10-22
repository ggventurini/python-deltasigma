python-deltasigma
===============

The **python-deltasigma** project aims to provide a 1:1 Python replacement of Richard 
Schreier's *excellent* **[MATLAB Delta Sigma Toolbox](http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox)**, upon which it is very heavily based. The original toolbox is also known as *delsig* or *delsigma*.

This project is at a very early stage. Take a look to files.ods for the current status. 

Pull requests welcome!

## Dependencies

Using **python-deltasigma** requires:

#### 1. **numpy** and **scipy** 

**numpy** and **scipy** can be installed through the steps described below.

On a Debian Linux system:

```
 # aptitude install python-numpy python-scipy python-matplotlib
```

Or, on a generic platform:

```
# pip install numpy scipy
```

More information on the [scipy install page](http://www.scipy.org/install.html).

#### 2. Slycot

**Slycot** is a *Python wrapper for selected SLICOT routines, notably including solvers for Riccati, Lyapunov and Sylvester equations.* (quoted from the project homepage.)

It is a dependency for **python-control** (see below). 

**Slycot** can be found at [repagh's Github repository](https://github.com/repagh/Slycot) and requires **numpy**, a **fortran compiler** such as gfortran and **BLAS/LAPACK 
libraries**. As, the README states, on a Debian Linux system, all of the above can be installed with:

```
# apt-get build-dep python-scipy
```

then check-out **Slycot** and install with distutils.

```
python setup.py install --user
```

If you are not using a Debian-based system, please check on the project page the dependencies to be installed.

#### 3. python-control

**python-control** can be installed downloading a release from [its homepage](http://sourceforge.net/projects/python-control/) or checking out its SVN repository with:

```
svn checkout svn://svn.code.sf.net/p/python-control/code/trunk python-control
```

Installing is straightforward with the distutils setup.py file:

```
$ python setup.py install --user
```

#### 4. matplotlib

**[matplotlib](http://matplotlib.org/)** is used for plotting and it is also very useful for visually inspecting your data.


#### 5. Extras

Building the documentation requires the **[sphinx](http://sphinx-doc.org/)** package.

If you plan to run the provided unit tests, then you should install **[setuptools](https://pypi.python.org/pypi/setuptools)**, used to access the reference function outputs. Testing *can* be automated with **[nose](https://pypi.python.org/pypi/nose/)**, issuing `$ nosetests -v pydelsigma/*.py`.

## Licensing and copyright notice

All original MATLAB code is Copyright (c) 2009, Richard Schreier. See the LICENSE file for the licensing terms.

The Python here provided is a derivative work from the above toolkit and subject to the same license terms.

The Python code is Copyright 2013, Giuseppe Venturini and the python-deltasigma contributors.
