# -*- coding: utf-8 -*-
# _padb.py
# This module provides the padb function.
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

"""This module provides the padb() function, which pads a matrix on the 
bottom.
"""

import numpy as np

def padb(x, n, val=0.):
    """Pad a matrix ``x`` on the bottom to length ``n`` with value ``val``.

    **Parameters:**

    x : ndarray
        The matrix to be padded.

    n : int
        The number of rows of the matrix after padding.

    val : scalar, optional
        The value to be used used for padding.

    .. note:: A 1-d array, for example ``a.shape == (N,)`` is reshaped to be
              a 1 column array: ``a.reshape((N, 1))``

    The empty matrix is assumed to be have 1 empty column.

    **Returns:**

    xp : 2-d ndarray
        The padded matrix.
    """
    if len(x.shape) == 1:
        xp = x.reshape((x.shape[0], 1))
    else:
        xp = x
    y = np.concatenate((xp, 
                        val*np.ones((n - xp.shape[0], xp.shape[1]))
                       ), axis=0)
    return y

