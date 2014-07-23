# -*- coding: utf-8 -*-
# _cancelPZ.py
# Module providing the cancelPZ function
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

"""Module providing the cancelPZ() function
"""

from __future__ import division
import numpy as np
import copy

from scipy.signal import lti

def cancelPZ(arg1, tol=1e-6):
    """Cancel zeros/poles in a SISO transfer function.

    **Parameters:**

    arg1 : arguments
        The `cancelPZ` function can be called with either 1, 2, 3 or 4 arguments.

    If one argument is used, it is a scipy `lti` object.

    If more arguments are used, they should be arranged in a tuple, the 
    following gives the number of elements in the tuple and their
    interpretation:

    * 2: (numerator, denominator)
    * 3: (zeros, poles, gain)
    * 4: (A, B, C, D)

    Each argument can be an array or sequence.

    tol : float, optional
        the absolute tolerance for pole, zero cancellation. Defaults to 1e-6.

    **Returns:**

    (z, p, k) : tuple
        A tuple containing zeros, poles and gain (unchanged) after poles, zeros 
        cancellation.

    """
    if not isinstance(arg1, lti):
        arg1 = lti(*arg1)
    z = copy.copy(arg1.zeros)
    p = copy.copy(arg1.poles)
    k = arg1.gain
    for i in range(max(z.shape) - 1, 0, -1):
        d = z[i] - p
        cancel = np.nonzero(np.abs(d) < tol)[0]
        if cancel.size:
            p = np.delete(p, cancel[0])
            z = np.delete(z, i)
    return z, p, k

