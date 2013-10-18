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

__author__ = "Giuseppe Venturini and the python-deltasigma contributors"
__copyright__ = "Copyright 2013, Giuseppe Venturini"
__credits__ = ["Giuseppe Venturini"]
__license__ = "BSD 3-Clause License"
__version__ = '0.0001alpha'
__maintainer__ = "Giuseppe Venturini"
__email__ = "ggventurini+github@gmail.com"
__status__ = "Pre-Pre-Alpha"

from ._SIunits import SIunits
from ._axisLabels import axisLabels
from ._bplogsmooth import bplogsmooth
from ._bquantize import bquantize
from ._constants import *
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
from ._impL1 import impL1
from ._infnorm import infnorm
from ._nabsH import nabsH
from ._padb import padb
from ._padl import padl
from ._padr import padr
from ._padt import padt
from ._plotPZ import plotPZ
from ._pulse import pulse
from ._rms import rms
from ._sinc_decimate import sinc_decimate
from ._thermometer import thermometer
from ._undbm import undbm
from ._undbp import undbp
from ._undbv import undbv
from ._utils import cplxpair, mfloor, rat, gcd, lcm, db
from ._zinc import zinc
