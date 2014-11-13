# -*- coding: utf-8 -*-
# _simulateDSM_cblas.py
# Module providing a cython implementation of simulateDSM,
# based on the BLAS library.
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
# The (future?) modifications are:
# Copyright (c) 2014, G. Venturini and the python-deltasigma contributors
#


# distutils: libraries = cblas

"""
Fast simulator for a generic delta sigma modulator using external cblas
=======================================================================
"""

import numpy as np
cimport numpy as np
import scipy as sp
__import__('scipy.signal')
from libc.math cimport floor, fabs

cdef extern from "cblas.h":
    enum CBLAS_ORDER:     CblasRowMajor, CblasColMajor
    enum CBLAS_TRANSPOSE: CblasNoTrans, CblasTrans, CblasConjTrans
    void cblas_dgemv(CBLAS_ORDER order, \
        CBLAS_TRANSPOSE TransA, int M, int N,\
        double alpha, double *A, int lda,\
        double *X, int incX,\
        double beta, double *Y, int incY)
    void cblas_dcopy(int N, double *X, int incX,\
        double *Y, int incY)

include '_simulateDSM_helper.pxi'

def simulateDSM(np.ndarray u, arg2, nlev=2, x0=0,
                int store_xn=False, int store_xmax=False, int store_y=False):

    # Make sure that nlev is a 1D int array
    cdef np.ndarray c_nlev
    try:
        c_nlev = np.asarray(nlev, dtype=np.int)
        if c_nlev.ndim > 1:
            raise TypeError()
        c_nlev=c_nlev.reshape(-1)
    except (ValueError, TypeError):
         raise ValueError(\
            "invalid argument: nlev must be convertible into a 1D int array")

    # Make sure that input is a matrix
    cdef np.ndarray c_u
    try:
        c_u = np.asarray(u, dtype=np.float64, order='C')
        if c_u.ndim > 2:
            raise TypeError()
        if c_u.ndim == 1:
            c_u = c_u.reshape(1, -1)
    except (ValueError, TypeError):
        raise ValueError(\
            "invalid argument: u must be convertible into a 2D float array")

    cdef int nu = c_u.shape[0]
    cdef int nq = c_nlev.shape[0]

    cdef int order

    try:
        if type(arg2) in (tuple, list) and len(arg2)==3:
            # Assume ntf in zpk form
            ntf_z=np.asarray(arg2[0], dtype=np.complex128)
            ntf_p=np.asarray(arg2[1], dtype=np.complex128)
            ntf_k=float(arg2[2])
            if ntf_z.ndim !=1 or ntf_p.ndim != 1:
                raise TypeError()
            form = 2
            order = ntf_z.shape[0]
        else:
            # Assume ABCD form
            ABCD = np.asarray(arg2, dtype=np.float64)
            if ABCD.ndim!=2:
                raise TypeError()
            if ABCD.shape[1] != nu+ABCD.shape[0]:
                raise TypeError()
            form = 1
            order = ABCD.shape[0]-nq
    except (ValueError, TypeError):
        raise ValueError('incorrect modulator specification')

    cdef np.ndarray c_x0, c_x0_temp
    # Assure that the state is a column vector
    try:
        if np.isscalar(x0) and x0 == 0:
            c_x0 = np.zeros((order, 1), dtype=np.float64)
        else:
            c_x0 = np.array(x0, dtype=np.float64, order='C')
            if c_x0.ndim < 1 or c_x0.ndim > 2:
                raise TypeError()
            c_x0=c_x0.reshape(-1, 1)
            if c_x0.shape[0]!=order:
                raise TypeError()
    except (ValueError, TypeError):
        raise ValueError('incorrect initial condition specification')
    c_x0_temp = np.empty_like(c_x0)

    cdef np.ndarray A, B1, B2, C, D1
    # Build ISO Model
    # note that B=hstack((B1, B2))
    if form == 1:
        A = np.asarray(ABCD[0:order, 0:order], dtype=np.float64, order='C')
        B1 = np.asarray(ABCD[0:order, order:order+nu],\
            dtype=np.float64, order='C')
        B2 = np.asarray(ABCD[0:order, order+nu:order+nu+nq],\
            dtype=np.float64, order='C')
        C = np.asarray(ABCD[order:order+nq, 0:order],\
            dtype=np.float64, order='C')
        D1 = np.asarray(ABCD[order:order+nq, order:order+nu], \
            dtype=np.float64, order='C')
    else:
        # Seek a realization of -1/H
        A, B2, C, D2 = sp.signal.zpk2ss(ntf_p, ntf_z, -1)
        C=C.real
        # Transform the realization so that C = [1 0 0 ...]
        Sinv = sp.linalg.orth(np.hstack((np.transpose(C), np.eye(order))))/ \
            np.linalg.norm(C)
        S = sp.linalg.inv(Sinv)
        C = np.dot(C, Sinv)
        if C[0, 0] < 0:
            S = -S
            Sinv = -Sinv
        A = np.asarray(S.dot(A).dot(Sinv), dtype=np.float64, order='C')
        B2 = np.asarray(np.dot(S, B2), dtype=np.float64, order='C')
        C = np.asarray(np.hstack(([[1.]], np.zeros((1,order-1)))),\
            dtype=np.float64, order='C')
        # C=C*Sinv;
        # D2 = 0;
        # !!!! Assume stf=1
        B1 = -B2
        D1 = np.asarray(1., dtype=np.float64)
        #B = np.hstack((B1, B2))

    # N is number of input samples to deal with
    cdef int N = c_u.shape[1]
    # v is output vector
    cdef np.ndarray v = np.empty((nq, N), dtype=np.float64)
    cdef np.ndarray y
    if store_y:
        # Need to store the quantizer input
        y = np.empty((nq, N), dtype=np.float64)
    else:
        y = np.empty((0,0), dtype=np.float64)
    cdef np.ndarray xn
    if store_xn:
        # Need to store the state information
        xn = np.empty((order, N), dtype=np.float64)
    cdef np.ndarray xmax
    if store_xmax:
        # Need to keep track of the state maxima
        xmax = np.abs(c_x0)
    else:
        xmax = np.empty(0, dtype=np.float64)

    # y0 is output before the quantizer
    cdef np.ndarray y0 = np.empty(nq, dtype=np.float64)

    cdef int i
    for i in xrange(N):
        # Compute y0 = np.dot(C, c_x0) + np.dot(D1, u[:, i])
        cblas_dgemv(CblasRowMajor, CblasNoTrans, nq, order,\
            1.0, dbldata(C), order, \
            dbldata(c_x0), 1, \
            0.0, dbldata(y0), 1)
        cblas_dgemv(CblasRowMajor, CblasNoTrans, nq, nu,\
            1.0, dbldata(D1), nu, \
            dbldata(c_u)+i, N, \
            1.0, dbldata(y0), 1)
        if store_y:
            #y[:, i] = y0[:]
            cblas_dcopy(nq, dbldata(y0), 1,\
            dbldata(y)+i, N)
        ds_quantize(nq, dbldata(y0), 1, \
            intdata(c_nlev), 1, \
            dbldata(v)+i, N)
        # Compute c_x0 = np.dot(A, c_x0) +
        #   np.dot(B, np.vstack((u[:, i], v[:, i])))
        cblas_dgemv(CblasRowMajor, CblasNoTrans, order, order,\
            1.0, dbldata(A), order, \
            dbldata(c_x0), 1,\
            0.0, dbldata(c_x0_temp), 1)
        cblas_dgemv(CblasRowMajor, CblasNoTrans, order, nu,\
            1.0, dbldata(B1), nu, \
            dbldata(c_u)+i, N, \
            1.0, dbldata(c_x0_temp), 1)
        cblas_dgemv(CblasRowMajor, CblasNoTrans, order, nq,\
            1.0, dbldata(B2), nq, \
            dbldata(v)+i, N, \
            1.0, dbldata(c_x0_temp), 1)
        # c_x0[:,1] = c_x0_temp[:,1]
        cblas_dcopy(order, dbldata(c_x0_temp), 1,\
            dbldata(c_x0), 1)
        if store_xn:
            # Save the next state
            #xn[:, i] = c_x0
            cblas_dcopy(order, dbldata(c_x0), 1,\
            dbldata(xn)+i, N)
        if store_xmax:
            # Keep track of the state maxima
            # xmax = np.max((np.abs(x0), xmax), 0)
            track_vabsmax(order, dbldata(xmax), 1,\
                dbldata(c_x0), 1)
    if not store_xn:
        xn = c_x0
    return v.squeeze(), xn.squeeze(), xmax, y.squeeze()

