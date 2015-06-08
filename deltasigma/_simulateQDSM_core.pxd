# -*- coding: utf-8 -*-
# _simulateQDSM.py
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

"""Module providing the simulateQDSM() core function
"""

import numpy
cimport numpy as np
import cython

#@cython.wraparound(False)
#@cython.nonecheck(False)
#@cython.boundscheck(False)
@cython.locals(N=cython.int, k=cython.float, v=np.ndarray,
               y=np.ndarray,xn=np.ndarray,xmax=np.ndarray)
 
#def simulateQDSM_core(u, A, B, C, D1, order, nlev, nq, x0):
cpdef inline simulateQDSM_core(np.ndarray[complex, ndim=2] u,
                               np.ndarray[complex, ndim=2] A,
                               np.ndarray[complex, ndim=2] B,
                               np.ndarray[complex, ndim=2] C,
                               np.ndarray[complex, ndim=2] D1,
                               int order, nlev, int nq,
                               np.ndarray[complex, ndim=2] x0)

@cython.locals(v=np.ndarray, ytmp=np.ndarray)
cdef inline ds_qquantize(np.ndarray[complex, ndim=1] y, n)

