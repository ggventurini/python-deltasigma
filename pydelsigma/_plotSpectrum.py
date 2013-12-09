# -*- coding: utf-8 -*-
# _plotSpectrum.py
# Module providing the plotSpectrum function
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

"""Module providing the plotSpectrum() function
"""

import pylab as plt

from ._logsmooth import logsmooth

def plotSpectrum(X, fin, fmt='-'):
    """Plot a smoothed spectrum on the current figure.
    """
    f, p = logsmooth(X, fin)
    plt.semilogx(f, p, fmt)
