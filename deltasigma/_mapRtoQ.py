# -*- coding: utf-8 -*-
# _mapRtoQ.py
# Module providing the mapRtoQ function
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

"""Module providing the mapRtoQ() function
"""

from __future__ import division
import numpy as np

def mapRtoQ(ABCDr):
    """Map a real ABCD matrix to a quadrature one.

    **Parameters:**

    ABCDr : ndarray
        A real matrix describing a quadrature system.

    ``ABCDr`` has its states paired (real, imaginary).

    **Returns:**

    (ABCDq, ABCDp) : tuple

    Where:

    ABCDq : ndarray
        is the quadrature (complex) version of ABCDr.
    ABCDp : ndarray
        is the mirror-image system matrix.

    .. note:: ``ABCDp`` is zero if ``ABCDr`` has no quadrature errors.
    """
    ABCD11 = ABCDr[::2,   ::2]
    ABCD12 = ABCDr[::2,  1::2]
    ABCD21 = ABCDr[1::2,  ::2]
    ABCD22 = ABCDr[1::2, 1::2]

    ABCDq = 0.5*(ABCD11 + ABCD22) + 0.5j*(ABCD21 - ABCD12);
    ABCDp = 0.5*(ABCD11 - ABCD22) + 0.5j*(ABCD21 + ABCD12);

    return ABCDq, ABCDp

