# -*- coding: utf-8 -*-
# _scaleABCD.py
# Module providing the scaleABCD function
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

"""Module providing the scaleABCD() function
"""

from __future__ import division, print_function
import numpy as np
import numpy.random as npr

from ._partitionABCD import partitionABCD
from ._simulateDSM import simulateDSM

def scaleABCD(ABCD, nlev=2, f=0, xlim=1, ymax=None, umax=None, N_sim=1e5, N0=10):
    """Scale the loop filter of a general delta-sigma modulator for dynamic range.

    The ABCD matrix is scaled so that the state maxima are less than the
    specified limits (``xlim``). As a side effect, the maximum stable input is
    determined in the process.

    **Parameters:**

    ABCD : ndarray
        The state-space description of the loop filter.

    nlev : int, optional
        The number of levels in the quantizer.

    f : scalar
        The normalized frequency of the test sinusoid.

    xlim : scalar or ndarray
        A vector or scalar specifying the limit for each state variable.

    ymax : scalar, optional
        The stability threshold. Inputs that yield quantizer inputs above ymax
        are considered to be beyond the stable range of the modulator.
        If not provided, it will be set to :math:`n_{lev} + 5`

    umax : scalar, optional
        The maximum allowable input amplitude. ``umax`` is calculated if it
        is not supplied.

    **Returns:**

    ABCDs : ndarray
        The state-space description of the scaled loop filter.

    umax : scalar
        The maximum stable input amplitude. Input sinusoids with amplitudes
        below this value should not cause the modulator states to exceed their
        specified limits.

    S : ndarray
        The diagonal scaling matrix S.

    `S` is defined such that::

        ABCDs = [[S*A*Sinv, S*B], [C*Sinv, D]]
        xs = S*x

    Where the multiplications are *matrix multiplications*.

    """
    if ymax is None:
        ymax = nlev + 5
    order = ABCD.shape[0] - 1
    xlim = xlim*np.ones((order,)) if np.isscalar(xlim) else xlim.reshape((-1, ))
    if np.isreal(ABCD).all():
        quadrature = False
    else:
        quadrature = True
    npr.seed(0) # So that this function is repeatable
    # Envelope for smooth start-up 
    raised_cosine = 0.5*(1 - np.cos(np.pi/N0*np.arange(N0)))

    if umax is None:
        # Simulate the modulator with DC or sine wave inputs to detect its stable
        # input range. 
        # First get a rough estimate of umax.
        ulist = np.arange(0.1, 1.1, 0.1)*(nlev - 1)
        umax = nlev - 1
        N = 1000.0
        u0 = np.hstack((
                        np.exp(2j*np.pi*f*np.arange(- N0, 0))*raised_cosine, \
                        np.exp(2j*np.pi*f*np.arange(0, N))
                      )) \
              + 0.01*np.dot(np.array([[1, 1j]]), npr.randn(2, N + N0))
        if not quadrature:
            u0 = np.real(u0)
        for u in ulist:
            if not quadrature:
                v, x, xmax, y = simulateDSM(u*u0, ABCD, nlev)
            else:
                raise NotImplementedError("simulateQDSM has not been implemented yet.")
                #v, x, xmax, y = simulateQDSM(u*u0, ABCD, nlev)
            if np.max(np.abs(y)) > ymax:
                umax = u
                # umax is the smallest input found which causes 'instability'
                break
        if umax == ulist[0]:
            msg = 'Modulator is unstable even with an input amplitude of %.1f.'\
             % umax
            raise RuntimeError(msg)
            
    # More detailed simulation
    N = N_sim
    u0 = np.hstack((
                    np.exp(2j*np.pi*f*np.arange(-N0, 0))*raised_cosine, \
                    np.exp(2j*np.pi*f*np.arange(0, N))
                  )) \
            + 0.01*np.dot(np.array([[1, 1j]]), npr.randn(2, N + N0))
    if not quadrature:
        u0 = np.real(u0)
    maxima = np.zeros((1, order)) - 1
    ulist = np.linspace(0.7*umax, umax, 10)
    for u in ulist:
        if not quadrature:
            v, x, xmax, y = simulateDSM(u*u0, ABCD, nlev)
        else:
            raise NotImplementedError("simulateQDSM has not been implemented yet.")
            #v, x, xmax, y = simulateQDSM(u*u0, ABCD, nlev)
        if np.max(np.abs(y)) > ymax:
            break
        umax = u
        maxima = np.max(np.vstack((maxima, xmax.T)), axis=0, keepdims=True)
    scale = (maxima/xlim).reshape((-1))
    S = np.diag(1.0/scale)
    Sinv = np.diag(scale)
    A, B, C, D = partitionABCD(ABCD)
    ABCDs = np.vstack((
                       np.hstack((np.dot(np.dot(S, A), Sinv), np.dot(S, B))),
                       np.hstack((np.dot(C, Sinv), D))
                     ))
    return ABCDs, umax, S

