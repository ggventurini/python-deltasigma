# -*- coding: utf-8 -*-
# _mod1.py
# Module providing the mod1 function
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

"""Module providing the mod1() utility function
"""

import numpy as np

from ._calculateTF import calculateTF

def mod1():
    """A description of the first-order modulator.

    **Returns:**

    ABCD, NTF, STF : a tuple of (ndarray, lti, lti)
        The elements are the ABCD matrix (ndarray), the NTF (lti object), the STF (lti object).

    """
    A = np.array([[1.]])
    B = np.array([[1., -1.]])
    C = np.array([[1.]])
    D = np.array([[0., 0.]])
    ABCD = np.vstack((np.hstack((A, B)), np.hstack((C, D))))
    H, G = calculateTF(ABCD)
    return ABCD, H, G

