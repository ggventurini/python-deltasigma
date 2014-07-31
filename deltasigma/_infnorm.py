# -*- coding: utf-8 -*-
# _infnorm.py
# This module provides the infnorm function.
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

"""This module provides the infnorm() function, which finds the infinity
norm of a z-domain transfer function.
"""

from __future__ import division
from warnings import warn
import numpy as np
from scipy.optimize import fminbound
from ._nabsH import nabsH
from ._evalTF import evalTF

def infnorm(H):
    """Find the infinity norm of a z-domain transfer function.

    **Parameters:**

    H : object
        the LTI description of the DT system, which can be in one of the
        following forms:

        * an LTI object,
        * a zpk tuple,
        * a (num, den) tuple,
        * an ABCD matrix (internally converted to zpk representation),
        * a list-like containing the A, B, C, D matrices (also internally
          converted to zpk representation).

    **Returns:**

    Hinf : float
           The infinity norm of ``H``.
    fmax : float
           The frequency to which Hinf corresponds.
    """
    # Get a rough idea of the location of the maximum.
    N = 129
    w = np.linspace(0, 2*np.pi, num=N, endpoint=True)
    dw = 2*np.pi/(N-1)
    Hval = evalTF(H, np.exp(1j*w))
    Hinf = np.max(np.abs(Hval))
    wi = np.where(np.abs(Hval) == Hinf)[0]

    # Home in using the scipy "fminbound" function.
    # original MATLAB code:
    #   wmax = fminbnd(nabsH, w(wi)-dw, w(wi)+dw, options, H);
    wmax = fminbound(nabsH, w[wi]-dw, w[wi]+dw, args=(H,), \
                     xtol=1e-08, maxfun=5000, full_output=0)

    if wmax is None:
        warn('Hinf: scipy.optimize.fminbound() failed.'
             + ' The result returned may not be very accurate.')
        wmax = w[wi]

    Hinf = -nabsH(wmax, H);
    fmax = wmax/(2*np.pi);
    # in the original Toolbox, wmax is returned (seems to be never used though)
    # rep? We return fmax.
    return Hinf, fmax
