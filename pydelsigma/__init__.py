# -*- coding: utf-8 -*-
# __init__.py
# python-deltasigma module init file
# Copyright 2013 Giuseppe Venturini
# This file is part of python-deltasigma.
#
# python-deltasigma is a 1:1 Python replacement of Richard Schreier's 
# MATLAB delta sigma toolbox (aka "delsigma"), upon which it is heavily based.
# The delta sigma toolbox is (c) 2009, Richard Schreier.
#
# python-deltasigma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# LICENSE file for the licensing terms.

"""
python-deltasigma
=================

The **MATLAB Delta Sigma Toolbox** with **0% MATLAB** and **a *lot* more
Python**.

The **python-deltasigma** is a Python package to *synthesize, simulate,
scale and map to implementable structures* **delta sigma modulators**.

It aims to provide **a 1:1 Python port** of Richard Schreier's
***excellent*** `MATLAB Delta Sigma
Toolbox <http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>`__,
the *de facto* standard tool for high-level delta sigma simulation, upon
which it is very heavily based.\ |githalytics.com alpha|

Status
------

This project is a work in progress. Not all functionality has been
ported. Take a look at
`files.csv <https://github.com/ggventurini/python-deltasigma/blob/master/files.csv>`__
for the current status.

|Build Status| |Coverage Status|

To have an idea of the currently implemented functionality, take a look
at the following ipython notebooks:

-  `dsdemo1 <http://nbviewer.ipython.org/gist/ggventurini/8040189>`__,
   notebook port of the interactive ``dsdemo1.m``.
-  `dsdemo2 <http://nbviewer.ipython.org/gist/ggventurini/8044644>`__,
   notebook port of the interactive ``dsdemo2.m``.
-  `dsdemo4 <http://nbviewer.ipython.org/gist/ggventurini/8255785/dsdemo4.ipynb>`__,
   notebook port of ``dsdemo4.m``. `Audio
   file <https://gist.github.com/ggventurini/8255785/raw/8fb7d94236b917e6d557fb538d3f35a3144c038c/sax.wav.b64>`__.
-  `dsexample1 <http://nbviewer.ipython.org/7251113>`__, python
   version of ``dsexample1.m``.
-  `dsexample2 <http://nbviewer.ipython.org/8323435>`__, python
   version of ``dsexample2.m``.
-  `dsexample3 <http://nbviewer.ipython.org/8323046>`__, python
   version of ``dsexample3.m``.

If you have some examples you would like to share, `send me a
mail <http://tinymailto.com/5310>`__, and I will add them to the above
list.

Further functionality is expected to be ported according to `the
ROADMAP <https://github.com/ggventurini/python-deltasigma/blob/master/ROADMAP.md>`__.

Install
-------

Supported platforms
~~~~~~~~~~~~~~~~~~~

python-deltasigma runs on every platform and arch. supported by its
dependencies:

-  *Platforms*: Linux, Mac OS X, Windows.

-  *Archs*: x86, x86\_64 and armf (arm with floating point unit).

Dependencies
~~~~~~~~~~~~

Using python-deltasigma requires **Python 2** or **3**, **numpy**,
**scipy** (>= 0.11.0) and **matplotlib**.

They are packaged by virtually all the major *Linux distributions*.

On a Debian Linux system, you may install them issuing:

::

     # aptitude install python python-numpy python-scipy python-matplotlib

Refer to your system documentation for more information.

On *Windows*, I hear good things about:

-  `Enthought Canopy <https://www.enthought.com/store/>`__, a Python
   distribution that carries both free and commercial versions, and

-  `Anaconda <https://store.continuum.io/cshop/anaconda/>`__, which
   offers its full version for free.

I do not run Windows, so I can't really provide more info (sorry),
except that people tell me they manage to have a working setup.

*Mac OS X* is also supported by `Enthought
Canopy <https://www.enthought.com/store/>`__ and
`Anaconda <https://store.continuum.io/cshop/anaconda/>`__, which likely
provide the easiest and fastest solution to have a scientific Python
stack up and running.

More information can be found on the `scipy install
page <http://www.scipy.org/install.html>`__ and on the `matplotlib
homepage <http://matplotlib.org/>`__.

I wrote in a different context some directions to `compile numpy and
scipy
yourself <https://github.com/ahkab/ahkab/wiki/Install:-numpy-and-scipy>`__,
which also apply here. Be warned, it can easily get complicated.

Extras
~~~~~~

Building the documentation requires the
`sphinx <http://sphinx-doc.org/>`__ package. It is an optional step,
as the `the latest documentation is available
online <http://python-deltasigma.readthedocs.org/en/latest/>`__, without
need for you to build it.

If you plan to run the provided unit tests, then you should install
`setuptools <https://pypi.python.org/pypi/setuptools>`__, used to
access the reference function outputs. Testing *can* be automated with
`nose <https://pypi.python.org/pypi/nose/>`__, issuing::

    $ nosetests -v pydelsigma/*.py

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

How to contribute
-----------------

Pull requests are welcome!

There are only a few *guidelines*, which can be overridden every time it
is reasonable to do so:

-  Please try to follow ``PEP8``.

-  Try to keep the functions signature identical. Parameters with
   ``NaN`` default values have their default value replaced with
   ``None``.

-  If a function has a varible number of return values, its Python port
   should implement the maximum number of return values.

Support python-deltasigma
~~~~~~~~~~~~~~~~~~~~~~~~~

*I do not want your money.* I develop this software because I enjoy it
and because I use it myself.

If you wish to support the development of ``python-deltasigma`` or you
find the package useful or you otherwise wish to contribute monetarily,
***please donate to cancer research instead:***

-  `Association for International Cancer Research
   (eng) <http://www.aicr.org.uk/donate.aspx>`__, or
-  `Fond. IRCCS Istituto Nazionale dei Tumori
   (it) <http://www.istitutotumori.mi.it/modules.php?name=Content&pa=showpage&pid=24>`__.

Consider `sending me a mail <http://tinymailto.com/5310>`__ afterwards,
***it makes for great motivation!***

Why this project was born
~~~~~~~~~~~~~~~~~~~~~~~~~

I like challenges, Delta Sigma modulation and I don't have the money for
my own MATLAB license.

With this Python package you can simulate Delta Sigma modulators for
free, on any PC.

I hope you find it useful.

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

.. |githalytics.com alpha| image:: https://cruel-carlota.pagodabox.com/36f25accf60f391456efe66910bf84f8
   :target: http://githalytics.com/ggventurini/python-deltasigma
.. |Build Status| image:: https://travis-ci.org/ggventurini/python-deltasigma.png?branch=master
   :target: https://travis-ci.org/ggventurini/python-deltasigma
.. |Coverage Status| image:: https://coveralls.io/repos/ggventurini/python-deltasigma/badge.png?branch=master
   :target: https://coveralls.io/r/ggventurini/python-deltasigma?branch=master

Key Functions
-------------

.. autosummary::
    :nosignatures:

    synthesizeNTF
    clans
    synthesizeChebyshevNTF
    simulateDSM
    simulateSNR
    realizeNTF
    stuffABCD
    mapABCD
    scaleABCD
    calculateTF
    realizeNTF_ct
    mapCtoD
    evalTFP

Functions for quadrature systems
--------------------------------

.. autosummary::
    :nosignatures:

    mapQtoR
    mapRtoQ

Other selected functions
------------------------

.. autosummary::
    :nosignatures:

    mod1
    mod2
    calculateSNR
    predictSNR
    partitionABCD
    infnorm
    impL1
    l1norm
    pulse
    rmsGain

General utilities for data processing
-------------------------------------

.. autosummary::
    :nosignatures:

    db
    dbm
    dbp
    dbv
    undbm
    undbp
    undbv
    ds_hann
    rms
    padb
    padl
    padr
    padt
    cplxpair
    mfloor
    rat
    gcd
    lcm

Plotting and data display utilitites
------------------------------------

.. autosummary::
    :nosignatures:

    plotPZ
    plotSpectrum
    figureMagic
    circ_smooth
    logsmooth
    DocumentNTF
    PlotExampleSpectrum
    axisLabels
    bilogplot
    bplogsmooth
    lollipop
    changeFig
    pretty_lti
    SIunits


Utility functions for Delta Sigma simulation
--------------------------------------------

.. autosummary::
    :nosignatures:

    bquantize
    bunquantize
    cancelPZ
    circshift
    delay
    ds_f1f2
    ds_freq
    ds_optzeros
    ds_quantize
    ds_synNTFobj1
    dsclansNTF
    evalMixedTF
    evalRPoly
    evalTF
    frespF1
    nabsH
    peakSNR
    sinc_decimate
    zinc

Detailed documentation
----------------------

"""

__author__ = "Giuseppe Venturini and the python-deltasigma contributors"
__copyright__ = "Copyright 2013, Giuseppe Venturini"
__credits__ = ["Giuseppe Venturini"]
__license__ = "BSD 3-Clause License"
__version__ = '0.0001alpha'
__maintainer__ = "Giuseppe Venturini"
__email__ = "ggventurini+github@gmail.com"
__status__ = "Pre-Pre-Alpha"

# Package testing can be done remotely, without display. This would make
# matplotlib fail (and consequently, the test itself).
# We check for $DISPLAY, but this makes us probably lose in portability, 
# does Windows have this environment variable defined?
# Then again, who runs the test suit on a head-less windows machine...
# ... for the time being the following should be OK. If in the future that
# feature is needed by somebody, we can switch to
# if not os.system('python -c "import matplotlib.pyplot as plt;plt.figure()"')
import matplotlib, os
if not ('DISPLAY' in os.environ or os.environ.get('READTHEDOCS', None)):
    matplotlib.use('Agg')

from ._DocumentNTF import DocumentNTF
from ._PlotExampleSpectrum import PlotExampleSpectrum
from ._SIunits import SIunits
from ._axisLabels import axisLabels
from ._bilogplot import bilogplot
from ._bplogsmooth import bplogsmooth
from ._bquantize import bquantize
from ._bunquantize import bunquantize
from ._calculateSNR import calculateSNR
from ._calculateTF import calculateTF
from ._cancelPZ import cancelPZ
from ._changeFig import changeFig
from ._circ_smooth import circ_smooth
from ._clans import clans
from ._constants import eps
from ._db import db
from ._dbm import dbm
from ._dbp import dbp
from ._dbv import dbv
from ._delay import delay
from ._ds_f1f2 import ds_f1f2
from ._ds_freq import ds_freq
from ._ds_hann import ds_hann
from ._ds_optzeros import ds_optzeros
from ._ds_quantize import ds_quantize
from ._ds_synNTFobj1 import ds_synNTFobj1
from ._dsclansNTF import dsclansNTF
from ._evalMixedTF import evalMixedTF
from ._evalRPoly import evalRPoly
from ._evalF0 import evalF0
from ._evalF1 import evalF1
from ._evalTF import evalTF
from ._evalTFP import evalTFP
from ._figureMagic import figureMagic
from ._frespF1 import frespF1
from ._impL1 import impL1
from ._infnorm import infnorm
from ._l1norm import l1norm
from ._logsmooth import logsmooth
from ._lollipop import lollipop
from ._mapABCD import mapABCD
from ._mapCtoD import mapCtoD
from ._mapQtoR import mapQtoR
from ._mapRtoQ import mapRtoQ
from ._mod1 import mod1
from ._mod2 import mod2
from ._nabsH import nabsH
from ._padb import padb
from ._padl import padl
from ._padr import padr
from ._padt import padt
from ._partitionABCD import partitionABCD
from ._peakSNR import peakSNR
from ._plotPZ import plotPZ
from ._plotSpectrum import plotSpectrum
from ._predictSNR import predictSNR
from ._pulse import pulse
from ._realizeNTF import realizeNTF
from ._realizeNTF_ct import realizeNTF_ct
from ._rms import rms
from ._rmsGain import rmsGain
from ._scaleABCD import scaleABCD
from ._simulateDSM import simulateDSM
from ._simulateSNR import simulateSNR
from ._sinc_decimate import sinc_decimate
from ._stuffABCD import stuffABCD
from ._synthesizeChebyshevNTF import synthesizeChebyshevNTF
from ._synthesizeNTF import synthesizeNTF
from ._thermometer import thermometer
from ._undbm import undbm
from ._undbp import undbp
from ._undbv import undbv
from ._utils import circshift, cplxpair, mfloor, pretty_lti, rat, gcd, lcm
from ._zinc import zinc
