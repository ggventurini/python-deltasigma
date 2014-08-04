# -*- coding: utf-8 -*-
# _rms.py
# This module provides the rms function.
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

"""This module provides the rms() function, which calculates the Root Mean 
Square (RMS) of a vector.
"""

import numpy as np
import numpy.linalg as la

def rms(x, no_dc=False):
    """Calculate the RMS value of ``x``.

    The Root Mean Square value of an array :math:`x` of length :math:`n` is defined as:

    .. math::

        x_{RMS} = \\sqrt{\\frac{1}{n}(x_1^2 + x_2^2 + ...+x_n^2)}

    **Parameters:**

    x : (N,) ndarray
        The input vector

    no_dc : boolean, optional
        If set to ``True``, the DC value gets subtracted from ``x`` first and the RMS is computed on the result.

    **Returns:**

    xrms : scalar
        as defined above

    """
    if no_dc:
        x = x - np.mean(x)
    return la.norm(x)/np.sqrt(max(x.shape))

