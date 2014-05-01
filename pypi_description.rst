It aims to provide **a 1:1 Python port** of Richard Schreier's
***excellent*** `MATLAB Delta Sigma
Toolbox <http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>`__,
the *de facto* standard tool for high-level delta sigma simulation, upon
which it is very heavily based.

Getting started
~~~~~~~~~~~~~~~

To have an idea of the currently implemented functionality, take a look
at these IPython notebooks: 

-  `dsdemo1 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsdemo1.ipynb>`__,
   notebook port of the interactive ``dsdemo1.m``.
-  `dsdemo2 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsdemo2.ipynb>`__,
   notebook port of the interactive ``dsdemo2.m``.
-  `dsdemo3 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsdemo3.ipynb>`__, 
   notebook port of the interactive `dsdemo3.m`.
-  `dsdemo4 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsdemo4.ipynb>`__,
   notebook port of ``dsdemo4.m``. `Audio
   file, right click to download <https://raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/sax.wav.b64>`__.
-  `dsexample1 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsexample1.ipynb>`__, python
   version of ``dsexample1.m``.
-  `dsexample2 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsexample2.ipynb>`__, python
   version of ``dsexample2.m``.
-  `dsexample3 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsexample3.ipynb>`__, python
   version of ``dsexample3.m``.

If you have some examples you would like to share, `send me a
mail <http://tinymailto.com/5310>`__, and I will add them to the above
list.

Further functionality is expected to be ported according to `the
ROADMAP <https://github.com/ggventurini/python-deltasigma/blob/master/ROADMAP.md>`__.

Documentation
-------------

1. You can find the included `package documentation
   online <http://python-deltasigma.readthedocs.org/en/latest/>`__.

2. The original MATLAB Toolbox provides in-depth documentation, which is
   very useful to understand what the toolbox is capable of. See
   `DSToolbox.pdf <https://github.com/ggventurini/python-deltasigma/blob/master/delsig/DSToolbox.pdf?raw=true>`__
   and
   `OnePageStory.pdf <https://github.com/ggventurini/python-deltasigma/blob/master/delsig/OnePageStory.pdf?raw=true>`__
   (*PDF warning*).

3. The book:

Richard Schreier, Gabor C. Temes, *Understanding Delta-Sigma Data
Converters*, ISBN: 978-0-471-46585-0, November 2004, Wiley-IEEE Press

is probably *the most authoritative resource on the topic*. Chapter 8-9
show how to use the MATLAB toolkit and the observations apply also to
this Python port. Links `on
amazon <http://www.amazon.com/Understanding-Delta-Sigma-Converters-Richard-Schreier/dp/0471465852>`__,
`on the Wiley-IEEE
press <http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0471465852,miniSiteCd-IEEE2.html>`__.

*I am not affiliated with neither the sellers nor the authors.*

Licensing and copyright notice
------------------------------

All original MATLAB code is Copyright (c) 2009, Richard Schreier. See
the LICENSE file for the licensing terms.

The Python code here provided is a derivative work from the above
toolkit and subject to the same license terms.

This package contains some source code from ``pydsm``, also based on the
same MATLAB toolbox. The ``pydsm`` package is copyright (c) 2012, Sergio
Callegari.

When not otherwise specified, the Python code is Copyright 2013,
Giuseppe Venturini and the python-deltasigma contributors.

MATLAB is a registered trademark of The MathWorks, Inc.
