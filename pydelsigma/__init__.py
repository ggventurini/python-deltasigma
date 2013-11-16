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

The **python-deltasigma** project aims to provide a 1:1 Python
replacement of Richard Schreier's MATLAB `delta sigma
toolbox <http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>`_.

**python-deltasigma** requires:
------------

Note: using virtualenv could be very useful to have an up-to-date setup without 
touching the packages installed by your OS.

1. **numpy** and **scipy**
^^^^^^^^^^^^^^^^^^^^^^^^^^

On a generic platform, they can be installed with:

::

    # pip install numpy scipy

More information on the `scipy install
page <http://www.scipy.org/install.html>`_.

2. Slycot
^^^^^^^^^

**Slycot** is a *Python wrapper for selected SLICOT routines, notably
including solvers for Riccati, Lyapunov and Sylvester equations.*
(quoted from the project homepage.)

It is a dependency for **python-control** (see below).

**Slycot** can be found at `repagh's Github
repository <https://github.com/repagh/Slycot>`_ and requires **numpy**,
a **fortran compiler** such as gfortran and **BLAS/LAPACK libraries**.
As, the README states, on a Debian Linux system, all of the above can be
installed with:

::

    # apt-get build-dep python-scipy

then check-out **Slycot** and install with distutils.

::

    python setup.py install --user

If you are not using a Debian-based system, please check on the project
page the dependencies to be installed.

3. python-control
^^^^^^^^^^^^^^^^^

**python-control** can be installed downloading a release from `its
homepage <http://sourceforge.net/projects/python-control/>`_ or checking
out its SVN repository with:

::

    svn checkout svn://svn.code.sf.net/p/python-control/code/trunk python-control

Installing is straightforward with the distutils setup.py file:

::

    $ python setup.py install --user

4. matplotlib
^^^^^^^^^^^^^

**`matplotlib <http://matplotlib.org/>`_** is used for plotting and it
is also very useful for visually inspecting your data.

5. Extras
^^^^^^^^^

Building the documentation requires the
**`sphinx <http://sphinx-doc.org/>`_** package.

If you plan to run the provided unit tests, then you should install
**`setuptools <https://pypi.python.org/pypi/setuptools>`_**, used to
access the reference function outputs. Testing *can* be automated with
**`nose <https://pypi.python.org/pypi/nose/>`_**, issuing
``$ nosetests -v pydelsigma/*.py``.

Licensing and copyright notice
------------------------------

All original MATLAB code is Copyright (c) 2009, Richard Schreier. See
the LICENSE file for the licensing terms.

The Python here provided is a derivative work from the above toolkit and
subject to the same license terms.

The Python code is Copyright 2013, Giuseppe Venturini and the
python-deltasigma contributors.
"""

__author__ = "Giuseppe Venturini and the python-deltasigma contributors"
__copyright__ = "Copyright 2013, Giuseppe Venturini"
__credits__ = ["Giuseppe Venturini"]
__license__ = "BSD 3-Clause License"
__version__ = '0.0001alpha'
__maintainer__ = "Giuseppe Venturini"
__email__ = "ggventurini+github@gmail.com"
__status__ = "Pre-Pre-Alpha"

from ._DocumentNTF import DocumentNTF
from ._PlotExampleSpectrum import PlotExampleSpectrum
from ._SIunits import SIunits
from ._axisLabels import axisLabels
from ._bplogsmooth import bplogsmooth
from ._bquantize import bquantize
from ._bunquantize import bunquantize
from ._calculateSNR import calculateSNR
from ._calculateTF import calculateTF
from ._circ_smooth import circ_smooth
from ._constants import eps
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
from ._dsclansObj import dsclansObj
from ._evalRPoly import evalRPoly
from ._evalTF import evalTF
from ._figureMagic import figureMagic
from ._impL1 import impL1
from ._infnorm import infnorm
from ._l1norm import l1norm
from ._lollipop import lollipop
from ._mapABCD import mapABCD
from ._mapQtoR import mapQtoR
from ._mapRtoQ import mapRtoQ
from ._nabsH import nabsH
from ._padb import padb
from ._padl import padl
from ._padr import padr
from ._padt import padt
from ._partitionABCD import partitionABCD
from ._peakSNR import peakSNR
from ._plotPZ import plotPZ
from ._predictSNR import predictSNR
from ._pulse import pulse
from ._realizeNTF import realizeNTF
from ._rms import rms
from ._rmsGain import rmsGain
from ._scaleABCD import scaleABCD
from ._simulateDSM import simulateDSM
from ._simulateSNR import simulateSNR
from ._sinc_decimate import sinc_decimate
from ._stuffABCD import stuffABCD
from ._synthesizeNTF import synthesizeNTF
from ._thermometer import thermometer
from ._undbm import undbm
from ._undbp import undbp
from ._undbv import undbv
from ._utils import cplxpair, mfloor, rat, gcd, lcm, db
from ._zinc import zinc
