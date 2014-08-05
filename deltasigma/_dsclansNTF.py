# -*- coding: utf-8 -*-
# _dsclansNTF.py
# Module providing the dsclansNTF function
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

"""Module providing the dsclansNTF() function
"""

import numpy as np
from ._utils import carray

def dsclansNTF(x, order, rmax, Hz):
    """ Conversion of clans parameters into a NTF.

    Translate x into H.
    I've changed the relationships between (zeta, wn) and x
    in order to guarantee LHP roots of the s-polynomial.
    
    Returns the NTF, a zpk tuple.
    """
    x = x.squeeze()
    Hz = carray(Hz)
    Hz = Hz.reshape((-1,))
    Hp = np.zeros((1,), dtype=np.complex128)
    odd = (order % 2 == 1)
    if odd:
        s = -x[0]**2.
        Hp[0] = rmax*(1. + s)/(1. - s)

    for i in range(0+1*odd, order, 2):
        Hp = np.hstack((Hp, np.zeros((2,))))
        zeta = x[i]**2
        wn = x[i + 1]**2
        s = np.roots(np.array((1, 2*zeta*wn, wn**2)))
        Hp[i:i+2] = rmax*(1. + s)/(1. - s)

    H = (Hz, Hp, 1.)
    return H

