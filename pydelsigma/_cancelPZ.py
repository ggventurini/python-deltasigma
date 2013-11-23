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

def cancelPZ(zpk1, tol=1e-6):
    """Cancel zeros/poles in a zpk system.
    """
    z, p, k = copy.copy(zpk1)
    for i in range(max(z.shape) - 1, 0, -1):
        d = z[i] - p
        cancel = np.nonzero(np.abs(d) < tol)[0]
        if cancel.size:
            p = np.delete(p, cancel[0])
            z = np.delete(z, i)
    return z, p, k

def test_cancelPZ():
    """Test unit for cancelPZ"""
    zt = np.array((1, 2, 3))
    pt = np.array((1 + 2e-5, 2 + .5e-5, 3))
    kt = 2
    zpkt = (zt, pt, kt)
    zr, pr, kr = cancelPZ(zpkt, tol=1e-5)
    np.allclose(zr, (1., ), atol=1e-8, rtol=1e-6)
    np.allclose(pr, (1. + 2e-5, ), atol=1e-8, rtol=1e-6)
    np.allclose(kr, 2., atol=1e-8, rtol=1e-6)
