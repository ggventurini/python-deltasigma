python-deltasigma
=================

A port of the **MATLAB Delta Sigma Toolbox** based on free software and
very little sleep

The **python-deltasigma** is a Python package to *synthesize, simulate,
scale and map to implementable topologies* **delta sigma modulators**.

It aims to provide **a 1:1 Python port** of Richard Schreier's
***excellent*** **`MATLAB Delta Sigma
Toolbox <http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>`__**,
the *de facto* standard tool for high-level delta sigma simulation, upon
which it is very heavily based. |githalytics.com alpha|

--------------

**Homepage:** `python-deltasigma.io <http://python-deltasigma.io>`__

**Documentation:**
`docs.python-deltasigma.io <http://docs.python-deltasigma.io>`__

**Latest version:** `0.1 <https://pypi.python.org/pypi/deltasigma/>`__

|Build Status| |Coverage Status| |PyPi version| |PyPi downloads| |BSD 2
clause license| |DOI BADGE|

Install
-------

``python-deltasigma`` runs on Linux, Mac OS X and Windows.

Installing requires **Python 2.6+** or **3.3+**, **numpy**, **scipy**
(>= 0.11.0) and **matplotlib**.

Strongly recommended: **Cython** - for significantly faster delta sigma
modulator simulations.

They are packaged by virtually all the major Linux distributions.

I do not run Windows, so I can't really provide more info (sorry),
except that people tell me they manage to have a working setup.

When the dependencies are satisfied, run:

::

    pip install deltasigma

to install the latest stable version from the `Python Package Index
(PYPI) <http://pypi.python.org>`__, or:

::

    python setup.py install

if you're installing a development version from Git.

Extras
~~~~~~

Install the **`sphinx <http://sphinx-doc.org/>`__** package to build the
documentation yourself.

The test suite requires
**`setuptools <https://pypi.python.org/pypi/setuptools>`__**, used to
access the reference function outputs.

Testing can be automated with
**`nose <https://pypi.python.org/pypi/nose/>`__**, issuing:

::

    nosetests -v deltasigma/*.py

Documentation
-------------

In addition to the notebooks found in the ``examples/`` directory,
ported from the MATLAB Delta Sigma toolbox:

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

   is probably *the most authoritative resource on the topic*. Chapter
   8-9 show how to use the MATLAB toolkit and the observations apply
   also to this Python port. Links `on
   amazon <http://www.amazon.com/Understanding-Delta-Sigma-Converters-Richard-Schreier/dp/0471465852>`__,
   `on the Wiley-IEEE
   press <http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0471465852,miniSiteCd-IEEE2.html>`__.

   *I am not affiliated with neither the sellers nor the authors.*

Licensing and copyright notice
------------------------------

All original MATLAB code is Copyright (c) 2009, Richard Schreier. See
the LICENSE file for the licensing terms.

The Python code here provided is a derivative work from the above
toolkit - a faithful port - and subject to the same license terms.

Credit goes to Richard Schreier for writing the MATLAB Delta Sigma
Toolbox, devising algorithms and implementations and doing an excellent
work at it.

This package contains some source code from ``pydsm``, also based on the
same MATLAB toolbox. The ``pydsm`` package is copyright (c) 2012, Sergio
Callegari.

When not otherwise specified, the Python code is Copyright 2013,
Giuseppe Venturini and the python-deltasigma contributors.

MATLAB is a registered trademark of The MathWorks, Inc.

.. |githalytics.com alpha| image:: https://cruel-carlota.pagodabox.com/36f25accf60f391456efe66910bf84f8
   :target: http://githalytics.com/ggventurini/python-deltasigma
.. |Build Status| image:: https://travis-ci.org/ggventurini/python-deltasigma.png?branch=master
   :target: https://travis-ci.org/ggventurini/python-deltasigma
.. |Coverage Status| image:: https://coveralls.io/repos/ggventurini/python-deltasigma/badge.png?branch=master
   :target: https://coveralls.io/r/ggventurini/python-deltasigma?branch=master
.. |PyPi version| image:: http://img.shields.io/badge/version-0.1-brightgreen.png
   :target: https://pypi.python.org/pypi/deltasigma/
.. |PyPi downloads| image:: https://pypip.in/download/deltasigma/badge.png
   :target: https://pypi.python.org/pypi/deltasigma/
.. |BSD 2 clause license| image:: http://img.shields.io/badge/license-BSD-brightgreen.png
   :target: https://raw.githubusercontent.com/ggventurini/python-deltasigma/master/LICENSE
.. |DOI BADGE| image:: https://zenodo.org/badge/doi/10.5281/zenodo.11167.png
   :target: http://dx.doi.org/10.5281/zenodo.11167
