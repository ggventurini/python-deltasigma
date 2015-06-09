# -*- coding: utf-8 -*-
# _calculateQTF.py
# Module providing the calculateQTF function
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

"""Module providing the calculateQTF() function
"""

from __future__ import division

import numpy as np

from scipy.signal import ss2tf, tf2zpk

from ._cancelPZ import cancelPZ
from ._partitionABCD import partitionABCD


def calculateQTF(ABCDr):
    """Calculate noise and signal transfer functions for a quadrature modulator

    **Parameters:**

    ABCDr : ndarray
        The ABCD matrix, in real form. You may call :func:`mapQtoR` to convert
        an imaginary (quadrature) ABCD matrix to a real one.

    **Returns:**

    ntf, stf, intf, istf : tuple of zpk tuples
        The quadrature noise and signal transfer functions.

    :raises RuntimeError: if the supplied ABCD matrix results in denominator mismatches.
    """
    A, B, C, D = partitionABCD(ABCDr, 4)

    #Construct an ABCD description of the closed-loop system
    # sys is a tuple in A, B, C, D form
    Acl = A + np.dot(B[:, 2:4], C)
    Bcl = B
    Ccl = C
    Dcl = np.hstack((D[:, 0:2], np.eye(2)))
    #sys = (A + np.dot(B[:, 2:4], C), B, C,
    #       np.hstack((D[:, 0:2], np.eye(2))))

    #Calculate the 2x4 matrix of transfer functions
    tfs = np.empty((2, 4), dtype=object)
    # Each tf is a tuple in num, den form
    # tf[i, j] corresponds to the TF from input j to output i
    for i in range(2):
        for j in range(4):
            tfs[i, j] = ss2tf(Acl, Bcl, Ccl[i, :], Dcl[i, :], input=j)

    #Reduce these to NTF, STF, INTF and ISTF
    if any(tfs[0, 2][1] != tfs[1, 3][1]):
        raise RuntimeError('TF Denominator mismatch. Location 1')

    ntf_x = (0.5 * (tfs[0, 2][0] + tfs[1, 3][0]), tfs[0, 2][1])
    intf_x = (0.5 * (tfs[0, 2][0] - tfs[1, 3][0]), tfs[0, 2][1])

    if any(tfs[0, 3][1] != tfs[1, 2][1]):
        raise RuntimeError('TF Denominator mismatch. Location 2')

    ntf_y = (0.5 * (tfs[1, 2][0] - tfs[0, 3][0]), tfs[1, 2][1])
    intf_y = (0.5 * (tfs[1, 2][0] + tfs[0, 3][0]), tfs[1, 2][1])

    if any(ntf_x[1] != ntf_y[1]):
        raise RuntimeError('TF Denominator mismatch. Location 3')
    if any(tfs[0, 0][1] != tfs[1, 1][1]):
        raise RuntimeError('TF Denominator mismatch. Location 4')

    stf_x = (0.5 * (tfs[0, 0][0] + tfs[1, 1][0]), tfs[0, 0][1])
    istf_x = (0.5 * (tfs[0, 0][0] - tfs[1, 1][0]), tfs[0, 0][1])

    if any(tfs[0, 1][1] != tfs[1, 0][1]):
        raise RuntimeError('TF Denominator mismatch. Location 5')

    stf_y = (0.5 * (tfs[1, 0][0] - tfs[0, 1][0]), tfs[1, 0][1])
    istf_y = (0.5 * (tfs[1, 0][0] + tfs[0, 1][0]), tfs[1, 0][1])

    if any(stf_x[1] != stf_y[1]):
        raise RuntimeError('TF Denominator mismatch. Location 6')

    # suppress warnings about complex TFs
    #warning('off')
    ntf = cancelPZ(tf2zpk(ntf_x[0] + 1j*ntf_y[0], ntf_x[1]))
    intf = cancelPZ(tf2zpk(intf_x[0] + 1j*intf_y[0], intf_x[1]))
    stf = cancelPZ(tf2zpk(stf_x[0] + 1j*stf_y[0], ntf_x[1]))
    istf = cancelPZ(tf2zpk(istf_x[0] + 1j*istf_y[0], intf_x[1]))
    #warning('on')

    return ntf, stf, intf, istf
