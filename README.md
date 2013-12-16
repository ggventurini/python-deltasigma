python-deltasigma
===============

The MATLAB Delta Sigma Toolbox with 0% MATLAB and a lot more Python.

The **python-deltasigma** project aims to provide **a 1:1 Python replacement** of Richard 
Schreier's *excellent* **[MATLAB Delta Sigma Toolbox](http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox)**, upon which it is very heavily based.[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/36f25accf60f391456efe66910bf84f8 "githalytics.com")](http://githalytics.com/ggventurini/python-deltasigma)

### Status

This project is at **an early stage**. Not all functionality has been ported and currently *most of the example files do not run*. Take a look at [files.csv](https://github.com/ggventurini/python-deltasigma/blob/master/files.csv) for the current status. 

**NEWS!** To have an idea of the currently implemented functionality, take a look at these **[preliminary results](http://nbviewer.ipython.org/7251113)**, which showcase `dsexample1.m`.

## Dependencies

Using **python-deltasigma** requires **numpy**, **scipy** (>= 0.11.0) and **matplotlib**.

They are packaged by virtually all the major Linux distributions. 

On a Debian Linux system, you may install them issueing:

```
 # aptitude install python-numpy python-scipy python-matplotlib
```

Refer to your system documentation for more information. 

On Windows, I hear good things about 
[Enthought Canopy](https://www.enthought.com/store/), a Python distribution 
that carries both free and commercial versions. I do not run Windows, so I 
can't really provide more info (sorry), except that people manage to have
a working setup. 

More information can be found on the 
[scipy install page](http://www.scipy.org/install.html) and on the 
[matplotlib homepage](http://matplotlib.org/).


#### Extras

Building the documentation requires the **[sphinx](http://sphinx-doc.org/)** package.

If you plan to run the provided unit tests, then you should install 
**[setuptools](https://pypi.python.org/pypi/setuptools)**, used to access the 
reference function outputs. Testing *can* be automated with 
**[nose](https://pypi.python.org/pypi/nose/)**, issuing 
`$ nosetests -v pydelsigma/*.py`.

## Documentation

1. You can find the [package documentation online](http://python-deltasigma.readthedocs.org/en/latest/pydelsigma.html).

2. The original MATLAB Toolbox provides good documentation, which is very useful to understand what the toolbox is capable of. See [DSToolbox.pdf](https://github.com/ggventurini/python-deltasigma/blob/master/delsig/DSToolbox.pdf?raw=true) and [OnePageStory.pdf](https://github.com/ggventurini/python-deltasigma/blob/master/delsig/OnePageStory.pdf?raw=true) (*PDF warning*).

3. The book:

Richard Schreier, Gabor C. Temes, *Understanding Delta-Sigma Data Converters*, ISBN: 978-0-471-46585-0, November 2004, Wiley-IEEE Press 

is an excellent resource on the topic. Links [on amazon](http://www.amazon.com/Understanding-Delta-Sigma-Converters-Richard-Schreier/dp/0471465852), [on the Wiley-IEEE press](http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0471465852,miniSiteCd-IEEE2.html). 

*I am not affiliated with neither the sellers nor the authors.*

## How to contribute

Pull requests are welcome!

There are only a few *guidelines*, which can be overridden every time it is reasonable to do so:

* Please try to follow `PEP8`. Except you are free to indent with tabs or spaces as you please (but please stick with your choice). 

* Try to keep the functions signature identical. Parameters with `NaN` default values have their default value replaced with `None`. 

* If a function has a varible number of return values, its Python port should implement the maximum number of return values.

## Licensing and copyright notice

All original MATLAB code is Copyright (c) 2009, Richard Schreier. 
See the LICENSE file for the licensing terms.

The Python code here provided is a derivative work from the above toolkit and 
subject to the same license terms.

This package contains some source code from `pydsm`, also based on the same 
MATLAB toolbox. The `pydsm` package is copyright (c) 2012, Sergio Callegari.

When not otherwise specified, the Python code is Copyright 2013, Giuseppe 
Venturini and the python-deltasigma contributors.

MATLAB is a registered trademark of The MathWorks, Inc.
