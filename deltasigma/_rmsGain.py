# -*- coding: utf-8 -*-
# _rmsGain.py
# Module providing the rmsGain function
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

"""Module providing the rmsGain() function
"""

import numpy as np
from scipy.linalg import norm
from ._evalTF import evalTF

def rmsGain(H, f1, f2, N=100):
    """Compute the root mean-square gain of a discrete-time TF.

    The computation is carried out over the frequency band ``(f1, f2)``,
    employing ``N`` discretization steps.

    **Parameters:**

    H : object
        The discrete-time transfer function. See  :func:`evalTF` for the supported types.

    f1 : scalar
        The start value in Hertz of the frequency band over which the gain is evaluated.

    f2 : scalar
        The end value (inclusive) in Hertz of the aforementioned frequency band.

    N : integer, optional
        The number of discretization points to be taken over specified interval.

    **Returns:**

    Grms : scalar
        The root mean-square gain
    """

    w = np.linspace(2*np.pi*f1, 2*np.pi*f2, N)
    g = norm(evalTF(H, np.exp(1j*w))) / np.sqrt(N)

    return g

