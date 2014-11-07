# -*- coding: utf-8 -*-
# _simulateDSM_helper.pyxi
# Helper inline functions for the cython implementations of simulateDSM
#
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
#
# This file is originally from `pydsm`. Little modifications have been
# performed - most prominently a bug in the ABCD matrix handling has been
# fixed. Many thanks to the original author.
#
# The original file is:
# Copyright (c) 2012, Sergio Callegari
# All rights reserved.
#
# The modifications are:
# Copyright (c) 2014, G. Venturini and the python-deltasigma contributors
#

cdef inline double dbl_sat(double x, double a, double b):
    return a if x <= a else b if x >= b else x

cdef inline void ds_quantize(int N, double* y, int y_stride, \
                             int* n, int n_stride, \
                             double* v, int v_stride):
    """Quantize a signal according to a given number of levels."""
    cdef int qi
    cdef double L
    for qi in range(N):
        if n[2*qi*n_stride] % 2 == 0:
            v[qi*v_stride] = 2*floor(0.5*y[qi*y_stride]) + 1
        else:
            v[qi*v_stride] = 2*floor(0.5*(y[qi*y_stride] + 1))
        L = n[2*qi*n_stride] - 1
        v[qi*v_stride] = dbl_sat(v[qi*v_stride], -L, L)

cdef inline void track_vabsmax(int N,\
                               double* vabsmax, int vabsmax_stride,\
                               double* x, int x_stride):
    cdef int i
    cdef double absx
    for i in range(N):
        absx = fabs(x[i*x_stride])
        if absx > vabsmax[i*vabsmax_stride]:
            vabsmax[i*vabsmax_stride] = absx

cdef inline double *dbldata(np.ndarray arr):
    return <double *>np.PyArray_DATA(arr)

cdef inline int *intdata(np.ndarray arr):
    return <int *>np.PyArray_DATA(arr)
