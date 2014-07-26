# -*- coding: utf-8 -*-
# _ds_quantize.py
# This module provides the ds_quantize function.
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

"""This module provides the ds_quantize() function, used to quantize signals
according to a user-specified quantizer characteristic.
"""

import numpy as np


def ds_quantize(y, n=2):
    """Quantize ``y``

    Quantize a vector :math:`y` to:

     * an odd integer in :math:`[-n+1, n-1]`, if :math:`n` is even, or
     * an even integer in :math:`[-n+1, n-1]`, if :math:`n` is odd.

    The quantizer implementation details are repeated here from its
    documentation for the user's convenience:

       The quantizer is ideal, producing integer outputs centered
       about zero. Quantizers with an even number of levels are of
       the mid-rise type and produce outputs which are odd integers.
       Quantizers with an odd number of levels are of the mid-tread
       type and produce outputs which are even integers.

    This definition gives the same step height for both mid-rise
    and mid-tread quantizers.

    .. image:: ../doc/_static/quantizer_model.png
        :align: center
        :alt: Quantizer model


    **Parameters:**

    n : int or ndarray, optional
        The number of levels in the quantizer. If ``n`` is an integer,
        then all rows of y are fed to the same quantizer.
        If ``n`` is a column vector, each of its elements specifies
        how to quantize the rows of ``y``.

    **Returns:**

    v : ndarray
        The quantized vector.

    .. seealso::
        :func:`bquantize`, :func:`bunquantize`

    """
    assert (np.round(n, 0) == n).all()  # did we get an int or an array of int?
    if not isinstance(n, np.ndarray):
        n = n * np.ones(y.shape)  # we got an int
    else:
        assert len(n.shape) == 1 or 1 in n.shape
        n = (np.ones((max(n.shape), y.shape[1])).T * n).T

    i = (n % 2 == 0)
    v = np.zeros(y.shape)
    v[i] = 2 * np.floor(0.5 * y[i]) + 1     # mid-rise quantizer
    v[~i] = 2 * np.floor(0.5 * (y[~i] + 1)) # mid-tread quantizer

    # Limit the output
    L = n - 1
    for m in (-1, 1):
        i = (m * v > L)
        if i.any():
            v[i] = m * L[i]
    return v

