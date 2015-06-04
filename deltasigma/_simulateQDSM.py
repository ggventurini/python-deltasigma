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

"""Module providing the simulateQDSM() function
"""

from __future__ import division, print_function

import copy

import numpy as np

import scipy

from scipy.signal import freqresp

from ._ds_quantize import ds_quantize
from ._evalTF import evalTF
from ._partitionABCD import partitionABCD
from ._utils import carray, diagonal_indices, _is_zpk, _is_A_B_C_D, _is_num_den

def simulateQDSM(u, arg2, nlev=2, x0=None):
    """Simulate a quadrature Delta-Sigma modulator

    This function computes the output of a quadrature delta-sigma modulator with
    input :math:`u`, a structure described by :math:`ABCD`, an initial state
    :math:`x_0` (default all zeros) and a quantizer whose number of levels is
    specified by :math:`n_{lev}`.

    For multiple quantizers, make :math:`n_{lev}` a 1D vector, for complex
    quantization to a diamond lattice, multiply :math:`n_{lev}` by :math:`j`.

    In that case, the shapes of the input parameters are:
    ``u.shape = (nu, N)``, ``nlev.shape = (nqi,)``,
    ``ABCD.shape = (order+nq, order+nq+nu)``.

    Alternatively, the modulator may be described by a supported NTF,
    representation, in particular it is recommended to use a zpk object. In this
    case, the STF is assumed to be 1.

    **Parameters:**

    u : ndarray
        The input signal to the modulator.
    arg2 : ndarray or a supported LTI representation
        A description of the modulator to simulate.
        An ndarray instance is interpreted as an ABCD description. Equivalently,
        the ABCD matrix may be supplied in ``(A, B, C, D)`` tuple form. All
        other supported modulator specifications result in a conversion to a zpk
        representation.
    nlev : int or sequence-like, optional
        The number of levels in the quantizer. If set to a sequence, each of the
        elements is assumed to be the number of levels associated with a
        quantizer. Defaults to ``2``.
    x0 : float or sequence-like, optional
        The initial states of the modulator. If it is set to a float, all states
        are assumed to have the same value, ``x0``. If it is set to a
        sequence-like objct (list, tuple, 1D ndarray and similar), each entry is
        assumed to be the value of one of the modulator states, in ascending
        order. Defaults to ``0``.

    **Returns:**

    v : ndarray
        The quantizer output.
    xn : ndarray
        The modulator states.
    xmax : nedarray
        The maximum value that each state reached during simulation.
    y : ndarray
        The quantizer input (ie the modulator output).
    """

    nu = u.shape[0]
    nlev = np.atleast_1d(nlev)
    nq = max(nlev.shape)
    if isinstance(arg2, scipy.signal.lti):
        k = arg2.k
        zeros = np.asarray(arg2.z)
        poles = np.asarray(arg2.p)
        form = 2
        order = max(zeros.shape)
    elif _is_zpk(arg2):
        zeros, poles, k = copy.deepcopy(arg2)
        zeros = np.asarray(zeros)
        poles = np.asarray(poles)
        form = 2
        order = max(zeros.shape)
    elif isinstance(arg2, np.ndarray): # ABCD
        if arg2.shape[1] > 2 and arg2.shape[1] == nu + arg2.shape[0]:
            # ABCD dimesions OK
            form = 1
            ABCD = arg2
            order = ABCD.shape[0] - nq
        else:
            raise ValueError('The ABCD argument does not have proper ' +
                             'dimensions.')
    elif _is_A_B_C_D(arg2):
        ABCD = np.vstack((np.hstack((np.atleast_2d(arg2[0]),
                                     np.atleast_2d(arg2[1]))),
                          np.hstack((np.atleast_2d(arg2[2]),
                                     np.atleast_2d(arg2[3])))))
        form = 1
        order = ABCD.shape[0] - nq
    elif _is_num_den(arg2):
        zeros, poles, k = tf2zpk(*arg)
        form = 2
        order = max(zeros.shape)
    else:
        raise TypeError('The second argument is neither an ABCD matrix nor ' +
                        'an NTF.')

    if x0 is None:
        x0 = np.zeros(shape=(order, ), dtype='float64')
    else:
        x0 = carray(x0)
    if form == 1:
        A, B, C, D = partitionABCD(ABCD, nq + nu)
        D1 = D[:, :nu]
    else:
        # Create a FF realization of 1-1/H.
        # Note that MATLAB's zp2ss and canon functions don't work for complex
        # TFs.
        A = np.zeros(shape=(order, order), dtype='float64')
        B2 = np.vstack((np.atleast_2d(1), np.zeros(shape=(order-1, 1),
                                                   dtype='float64')))
        diag = diagonal_indices(A, 0)
        A[diag] = zeros
        subdiag = diagonal_indices(A, -1)
        A[subdiag] = 1.
        # Compute C st C*inv(zI-A)*B = 1-1/H(z);
        w = 2*np.pi*np.random.rand(2*order)
        desired = 1 - 1.0/evalTF((zeros, poles, k), np.exp(1j*w))
        # suppress warnings about complex TFs ???
        sys = scipy.signal.lti(A, B2, np.eye(order), np.zeros((order, 1)))
        # not clear why this is being reshaped
        sysresp = np.reshape(freqresp(sys, w)[1], order, max(w.shape))
        C = desired/sysresp
        # !!!! Assume stf=1
        B1 = -B2
        B = np.hstack((B1, B2))
        D1 = 1
    N = max(u.shape)
    v = np.zeros(shape=(nq, N), dtype='float64')
    y = np.zeros(shape=(nq, N), dtype='float64')
    # Need to store the state information
    xn = np.zeros(shape=(order, N), dtype='float64')
    # Need to keep track of the state maxima
    xmax = abs(x0.copy())
    for i in range(N):
        y[:, i] = np.dot(C, x0) + np.dot(D1, u[:, i])
        v[:, i] = ds_qquantize(y[:, i], nlev)
        x0 = np.dot(A, x0) + np.dot(B, np.vstack((u[:, i], v[:, i])))
        # Save the next state
        xn[:, i] = x0
        # Keep track of the state maxima
        xmax = np.max(np.abs(x0), xmax)
    return v, xn, xmax, y

def ds_qquantize(y, n=2):
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

