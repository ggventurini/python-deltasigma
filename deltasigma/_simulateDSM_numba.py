# -*- coding: utf-8 -*-
# _simulateDSM_numba.py
# Module providing the Numba simulateDSM function
# Copyright 2014 Giuseppe Venturini & Shayne Hodge
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

"""Module providing the simulateDSM() function implemented in Numba
"""

import collections

from warnings import warn

import numpy as np

from _utils import carray

from numba import double, int32, complex64
from numba.types import pyobject
from numba.decorators import jit, autojit

# Using object mode, is nopython mode faster?  Array slicing doesn't work in it
# Input data type assumptions
# u = 2D double array
# arg2 = float tuple?
# nlev = 1D int array
# x0 = 1D double array

# Temporarily rewriting to work only the ABCD form

def simulateDSM(u, arg2, nlev=2, x0=0):
    if isinstance(arg2, np.ndarray):
        nlev = carray(nlev)
        u = np.array(u) if not hasattr(u, 'ndim') else u
        if not max(u.shape) == np.prod(u.shape):
            warn("Multiple input delta sigma structures have had little testing.")
        if u.ndim == 1:
            u = u.reshape((1, -1))
        ABCD = np.asarray(arg2, dtype=np.float64)
        nq = 1 if np.isscalar(nlev) else nlev.shape[0]
        if ABCD.shape[1] != ABCD.shape[0] + u.shape[0]:
            raise ValueError('The ABCD argument does not have proper dimensions.')
        if not isinstance(x0, collections.Iterable):
            x0 = x0*np.ones((ABCD.shape[0] - nq, 1))
        else:
            x0 = np.array(x0).reshape((-1, 1))
        print('u = {}'.format(u))
        print('ABCD = {}'.format(ABCD))
        print('nq = {}'.format(nq))
        print('nlev = {}'.format(nlev))
        print('x0 = {}'.format(x0))
        v, xn, xmax, y = simulateDSM_ABCD(u, ABCD, nq, nlev, x0)
        return v, xn, xmax, y
    else:
        raise ValueError("Stop trying to do things that aren't supported yet!")

@jit(argtypes=([double[:, :], double[:, :], int32, int32[:], double[:]]))
#@autojit
def simulateDSM_ABCD(u, ABCD, nq, nlev=2, x0=0):
    # Removing type checking that is conflicting with Numba
    #if (hasattr(arg2, 'inputs') and not arg2.inputs == 1) or \
    #   (hasattr(arg2, 'outputs') and not arg2.outputs == 1):
    #        raise TypeError("The supplied TF isn't a SISO transfer function.")
    nu = u.shape[0]
    order = ABCD.shape[0] - nq
    A = ABCD[:order, :order]
    B = ABCD[:order, order:order+nu+nq]
    C = ABCD[order:order+nq, :order]
    D1 = ABCD[order:order+nq, order:order+nu]
    N = u.shape[1]
    v = np.empty((nq, N), dtype=np.float64)
    y = np.empty((nq, N), dtype=np.float64)     # to store the quantizer input
    xn = np.empty((order, N), dtype=np.float64) # to store the state information
    xmax = np.abs(x0) # to keep track of the state maxima

    for i in range(N):
        # y0 needs to be cast to real because ds_quantize needs real
        # inputs. If quantization were defined for complex numbers,
        # this cast could be removed
        y0 = np.real(np.dot(C, x0) + np.dot(D1, u[:, i]))
        y[:, i] = y0
        v[:, i] = ds_quantize(y0, nlev)
        x0 = np.dot(A, x0) + np.dot(B, np.vstack((u[:, i], v[:, i])))
        xn[:, i] = np.real_if_close(x0.T)
        xmax = np.max(np.hstack((np.abs(x0), xmax)), axis=1, keepdims=True)

    return v.squeeze(), xn.squeeze(), xmax, y.squeeze()

#@jit(arg_types=[double[:, :], int32[:, :]])
@autojit
def ds_quantize(y, n):
    """v = ds_quantize(y,n)
    Quantize y to:

    * an odd integer in [-n+1, n-1], if n is even, or
    * an even integer in [-n, n], if n is odd.

    This definition gives the same step height for both mid-rise
    and mid-tread quantizers.
    """
    v = np.zeros(y.shape)
    for qi in range(n.shape[0]):
        if n[qi] % 2 == 0: # mid-rise quantizer
            v[qi, 0] = 2*np.floor(0.5*y[qi, 0]) + 1
        else: # mid-tread quantizer
            v[qi, 0] = 2*np.floor(0.5*(y[qi, 0] + 1))
        L = n[qi] - 1
        v[qi, 0] = np.sign(v[qi, 0])*np.min((np.abs(v[qi, 0]), L))
    return v

