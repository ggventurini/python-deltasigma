# -*- coding: utf-8 -*-
# _sinc_decimate.py
# This module provides the sinc_decimate function.
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

"""This module provides the sinc_decimate() function, which decimates a vector 
by a sinc filter of specified order and length.
"""

from __future__ import division
import numpy as np

def sinc_decimate(x, m, r):
    """Decimate ``x`` by an ``m``-th order sinc filter of length ``r``.
    """
    x = x[:]
    for _ in range(m):
        x = np.cumsum(x)
        x = np.concatenate((x[:r], x[r:] - x[:-r]), axis=0)/r
    return x[r-1::r]

