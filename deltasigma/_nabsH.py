# -*- coding: utf-8 -*-
# _nabsH.py
# This module provides the nabsH function.
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

"""This module provides the nabsH() function, which computes the negative of
the absolute value of H(z).
"""

import numpy as np

from ._evalTF import evalTF

def nabsH(w, H):
    """Computes the negative of the absolute value of H.

    The computation is performed at the specified angular
    frequency ``w``, on the unit circle.

    This function is used by :func:`infnorm`.
    """
    z = np.exp(1j*w)
    return -np.abs(evalTF(H, z))

