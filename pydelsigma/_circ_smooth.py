# -*- coding: utf-8 -*-
# _circ_smooth.py
# Module providing the circ_smooth function
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

"""Module providing the circ_smooth() function
"""

from __future__ import division
import numpy as np
from ._ds_hann import ds_hann
from ._utils import circshift

def circ_smooth(x, n):
    """function y = circ_smooth(x,n)"""
    assert len(x.shape) == 1 or 1 in x.shape
    assert n % 2 == 0
    nx = max(x.shape)
    w = ds_hann(n)/(n/2.)
    xw = np.convolve(x, w)
    yp = np.hstack((xw[n - 1:nx], xw[:n - 1] + xw[nx:]))
    y = circshift(yp,[int(n/2. - 1)])
    return y
    
def test_circ_smooth():
    """Test function for circ_smooth()
    """
    import pkg_resources
    from scipy.io import loadmat
    A = np.arange(1, 101)
    b = circ_smooth(A, 16)
    fname = pkg_resources.resource_filename(__name__, "test_data/test_circ_smooth.mat")
    bt = loadmat(fname)['b']
    assert np.allclose(bt, b, atol=1e-8, rtol=1e-5) 
