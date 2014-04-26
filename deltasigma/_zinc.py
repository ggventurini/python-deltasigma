# -*- coding: utf-8 -*-
# _zinc.py
# The zinc function.
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

"""This module provides the zinc() function which calculates the magnitude
response of a cascade of comb filters.
"""

import numpy as np


def zinc(f, m=64, n=1):
    """Calculate the magnitude response of a cascade of ``n`` ``m``-th order comb filters.

    The magnitude of the filter response is calculated mathematically as:

    .. math::

        \\left|H(f)\\right| = \\left|\\frac{\\mathrm{sinc}(m f)}{\\mathrm{sinc}(f)}\\right|^n

    **Parameters:**

    f : ndarray
        The frequencies f at which the magnitude response is evaluated.

    m : int, optional
        The order of the comb filters.

    n : int, optional
        The number of comb filters in the cascade.

    **Returns:**

    HM : ndarray
        The magnitude of the frequency response of the cascade filter.

    """
    return np.fabs(np.sinc(m * f) / np.sinc(f)) ** n


def test_zinc():
    """Test function for zinc()"""
    ref = [1.0000,    0.9985,    0.9941,    0.9867,    0.9765,    0.9635,
           0.9478,    0.9295,    0.9087,    0.8855,    0.8602,    0.8329,
           0.8038,    0.7730,    0.7408,    0.7074,    0.6729,    0.6377,
           0.6019,    0.5658,    0.5295,    0.4933,    0.4574,    0.4221,
           0.3874,    0.3536,    0.3208,    0.2892,    0.2590,    0.2302,
           0.2031,    0.1776,    0.1538,    0.1319,    0.1118,    0.0936,
           0.0772,    0.0626,    0.0499,    0.0389,    0.0295,    0.0217,
           0.0154,    0.0104,    0.0066,    0.0038,    0.0020,    0.0008,
           0.0002,    0.0000,    0.0000]
    f = np.arange(0, 0.51, 0.01)
    test = zinc(f, 2, 3)
    assert np.allclose(test, ref, atol=1e-4, rtol=1e-4)
