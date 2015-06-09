# -*- coding: utf-8 -*-
# _simulateQDSM_core.py
# Module providing the simulateQDSM function
# Copyright 2015 Giuseppe Venturini
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

"""Module providing the core of the simulateQDSM() function
"""

from __future__ import division, print_function

import numpy as np

from ._ds_quantize import ds_quantize

def simulateQDSM_core(u, A, B, C, D1, order, nlev, nq, x0):
    N = u.shape[1]
    v = np.zeros(shape=(nq, N), dtype='complex128')
    y = np.zeros(shape=(nq, N), dtype='complex128')
    # Need to store the state information
    xn = np.zeros(shape=(order, N), dtype='complex128')
    # Need to keep track of the state maxima
    xmax = abs(x0.copy())
    for i in range(N):
        y[:, i] = np.dot(C, x0) + np.dot(D1, u[:, i].reshape((-1, 1)))
        v[:, i] = ds_qquantize(y[:, i], nlev)
        x0 = np.dot(A, x0) + np.dot(B, np.vstack((u[:, i], v[:, i])))
        # Save the next state
        xn[:, i] = np.squeeze(x0)
        # Keep track of the state maxima
        xmax = np.max((np.abs(x0), xmax), axis=0)
    return v, xn, xmax, y

def ds_qquantize(y, n):
    """Quadrature quantization
    """
    if np.isreal(n):
        v = ds_quantize(np.real(y), n) + 1j*ds_quantize(np.imag(y), n)
    else:
        ytmp = np.concatenate((np.real(y) + np.imag(y),
                               np.real(y) - np.imag(y)))
        v = np.dot(ds_quantize(ytmp, np.abs(n)),
                   np.array([[1 + 1j], [1 - 1j]]))/2.
    return v

