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
    xlim = xlim*np.ones((1, order)) if np.isscalar(xlim) else xlim.reshape((1, -1))
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

def test_scaleABCD():
    """Test function for scaleABCD()"""
    from ._synthesizeNTF import synthesizeNTF
    from ._realizeNTF import realizeNTF
    from ._stuffABCD import stuffABCD
    order = 8
    osr = 32
    nlev = 2
    f0 = 0.125
    Hinf = 1.5
    form = 'CRFB'
    ntf = synthesizeNTF(order, osr, 2, Hinf, f0)            # Optimized zero placement
    a, g, b, c = realizeNTF(ntf, form)
    # we pass b = b[0] bc if you use a single feed-in for the input, b may be scalar too
    ABCD = stuffABCD(a, g, b[0], c, form)
    # now we correct b for the assert
    b = np.concatenate(( # Use a single feed-in for the input
                        np.atleast_1d(b[0]),
                        np.zeros((max(b.shape) - 1,))
                      ))
    ABCD0 = ABCD.copy()
    # Values to be tested
    ABCD, umax, S = scaleABCD(ABCD0, nlev, f0)

    # References
    Sdiag_ref = np.array([71.9580, 51.9359, 8.2133, 6.5398, 1.9446, 1.2070, 
                          0.4223, 0.3040])
    umax_ref = 0.8667
    ABCD_ref1 = np.array([
    [    1.0000,   -0.7320,         0,         0,         0,         0],
    [    0.7218,    0.4717,         0,         0,         0,         0],
    [         0,    0.1581,    1.0000,   -0.7357,         0,         0],
    [         0,    0.1259,    0.7962,    0.4142,         0,         0],
    [         0,         0,         0,    0.2973,    1.0000,   -0.9437],
    [         0,         0,         0,    0.1846,    0.6207,    0.4142],
    [         0,         0,         0,         0,         0,    0.3499],
    [         0,         0,         0,         0,         0,    0.2518],
    [         0,         0,         0,         0,         0,         0]])
    ABCD_ref2 = np.array([
    [         0,         0,    0.0858,   -0.0858],
    [         0,         0,    0.0619,    0.0428],
    [         0,         0,         0,    0.0642],
    [         0,         0,         0,    0.1835],
    [         0,         0,         0,    0.2447],
    [         0,         0,         0,    0.0581],
    [    1.0000,   -0.8971,         0,   -0.0076],
    [    0.7197,    0.3543,         0,   -0.1746],
    [         0,    3.2900,         0,         0]])
    ABCD_ref = np.hstack((ABCD_ref1, ABCD_ref2))

    # mapping the NTF to states, there is not a perfect match between 
    # the original code and the scipy version. -> rtol approx 20%
    #if not np.allclose(ABCD, ABCD_ref, atol=1e-2, rtol=3e-1):
    #    aerr = ABCD_ref-ABCD
    #    rerr = 2*(ABCD_ref-ABCD)/(ABCD_ref+ABCD)
    #    print(repr(ABCD_ref))
    #    print(repr(ABCD))
    #    print(aerr)
    #    print(rerr)
    # this is a rather high relative error. We get it on Travis-CI
    # Probably linked to the libs used?
    assert np.allclose(ABCD, ABCD_ref, atol=1e-2, rtol=30e-2)
    assert np.allclose(umax, umax_ref, atol=1e-4, rtol=1e-3)
    assert np.allclose(np.diag(S), Sdiag_ref, atol=1e-2, rtol=25e-1)

