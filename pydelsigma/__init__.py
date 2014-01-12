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
The **python-deltasigma** project aims to provide a 1:1 Python
replacement of Richard Schreier's MATLAB `delta sigma
toolbox <http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>`_.

Getting started
---------------

Dependencies
^^^^^^^^^^^^

This toolbox requires `numpy <http://www.numpy.org/>`_, 
`scipy <http://www.scipy.org/>`_ and 
`matplotlib <http://matplotlib.org/>`_ to be available on your system.

Extras
^^^^^^

Building the documentation requires the
`sphinx <http://sphinx-doc.org/>`_ package.

If you plan to run the provided unit tests, then you should install
`setuptools <https://pypi.python.org/pypi/setuptools>`_, used to
access the reference function outputs. 

Testing can be automated with
`nose <https://pypi.python.org/pypi/nose/>`_, issuing:

::

    $ nosetests -v pydelsigma/*.py

:note: using `virtualenv` could be very useful to have an up-to-date setup without affecting the packages installed by your OS.

Licensing and copyright notice
------------------------------

The MATLAB Delta-Sigma Toolbox and all original MATLAB code is 
copyright (c) 2009, Richard Schreier. See the LICENSE file for the licensing terms.

The Python code here provided is a derivative work from the above toolkit and
subject to the same license terms.

The Python code in this package is Copyright 2013, Giuseppe Venturini and the
python-deltasigma contributors, unless otherwise noted.

This package contains some source code from pydsm, a previous Python port of 
the MATLAB Delta Sigma Toolbox. Pydsm is copyright (c) 2012, Sergio Callegari.

Functions
---------

.. autosummary::
    :nosignatures:

    DocumentNTF
    PlotExampleSpectrum
    SIunits
    axisLabels
    bilogplot
    bplogsmooth
    bquantize
    bunquantize
    calculateSNR
    calculateTF
    cancelPZ
    changeFig
    circ_smooth
    circshift
    clans
    eps
    db
    dbm
    dbp
    dbv
    delay
    ds_f1f2
    ds_freq
    ds_hann
    ds_optzeros
    ds_quantize
    ds_synNTFobj1
    dsclansNTF
    evalMixedTF
    evalRPoly
    evalF0
    evalF1
    evalTF
    evalTFP
    figureMagic
    frespF1
    impL1
    infnorm
    l1norm
    logsmooth
    lollipop
    mapABCD
    mapCtoD
    mapQtoR
    mapRtoQ
    mod1
    mod2
    nabsH
    padb
    padl
    padr
    padt
    partitionABCD
    peakSNR
    plotPZ
    plotSpectrum
    predictSNR
    pulse
    realizeNTF
    realizeNTF_ct
    rms
    rmsGain
    scaleABCD
    simulateDSM
    simulateSNR
    sinc_decimate
    stuffABCD
    synthesizeChebyshevNTF
    synthesizeNTF
    thermometer
    undbm
    undbp
    undbv
    cplxpair
    mfloor
    rat
    gcd
    lcm
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
from ._utils import circshift, cplxpair, mfloor, rat, gcd, lcm
from ._zinc import zinc
