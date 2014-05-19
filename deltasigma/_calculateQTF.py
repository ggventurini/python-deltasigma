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

from ._partitionABCD import partitionABCD

def calculateQTF(ABCDr):
    """Calculate noise and signal transfer functions for a quadrature modulator

    **Parameters:**

    ABCDr : ndarray
        The ABCD matrix.

    **Returns:**

    [ntf, stf, intf, istf] : list of zpk tuples
        A list containing the quadrature noise and signal transfer functions
    """
    A, B, C, D = partitionABCD(ABCDr, 4)

    #Construct an ABCD description of the closed-loop system
    sys = (A + np.dot(B[:, 2:4], C), 
           B, 
           C, 
           np.array([D[:, 0:2], eye(2)]).reshape((1, -1)), 
           1)

    #Calculate the 2x4 matrix of transfer functions
    tfs = np.empty((2, 4), dtype=object)
    for i in range(2):
        for j in range(4):
            tfs[i, j] = (ss2tf(sys, input=i)[0][j], ss2tf(sys, input=i)[0])

    #Reduce these to NTF, STF, INTF and ISTF
    if any(tfs[0, 2][1] != tfs[1, 3][1]):
        error('TF Denominator mismatch. Location 1')

    ntf_x =  (0.5*(tfs[0, 2][0] + tfs[1, 3][0]), tfs[0, 2][1])
    intf_x = (0.5*(tfs[0, 2][0] - tfs[1, 3][0]), tfs[0, 2][1])

    if any(tfs[0, 3][1] != tfs[1, 2][1]):
        error('TF Denominator mismatch. Location 2')

    ntf_y = (0.5*(tfs[1, 2][0] - tfs[0, 3][0]), tfs[1, 2][1])
    intf_y = (0.5*(tfs[1, 2][0] + tfs[0, 3][0]), tfs[1, 2][1])

    if any(ntf_x[1] != ntf_y[1]):
        error('TF Denominator mismatch. Location 3')
    if any(tfs[0, 0][1] != tfs[1, 1][1]):
        error('TF Denominator mismatch. Location 4')

    stf_x =  (0.5*(tfs[0, 0][0] + tfs[1, 1][0]), tfs[0, 0][1])
    istf_x = (0.5*(tfs[0, 0][0] - tfs[1, 1][0]), tfs[0, 0][1])

    if any(tfs[0, 1][1] != tfs[1, 0][1]):
        error('TF Denominator mismatch. Location 5')

    stf_y = (0.5*(tfs[1, 0][0] - tfs[0, 1][0]), tfs[1, 0][1])
    istf_y = (0.5*(tfs[1, 0][0] + tfs[0, 1][0]), tfs[1, 0][1])

    if any(stf_x[1] != stf_y[1]):
        error('TF Denominator mismatch. Location 6')

    # suppress warnings about complex TFs
    #warning('off')
    ntf =  cancelPZ(tf2zpk(( ntf_x.num[0] + 1j* ntf_y.num[0],  ntf_x.den[0])))
    intf = cancelPZ(tf2zpk((intf_x.num[0] + 1j*intf_y.num[0], intf_x.den[0])))
    stf =  cancelPZ(tf2zpk(( stf_x.num[0] + 1j* stf_y.num[0],  ntf_x.den[0])))
    istf = cancelPZ(tf2zpk((istf_x.num[0] + 1j*istf_y.num[0], intf_x.den[0])))
    #warning('on')

    return ntf, stf, intf, istf
