# -*- coding: utf-8 -*-
# _bquantize.py
# Bipolar quantization module
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

"""This module provides the bquantize() function, used to bidirectionally 
quantize a vector to signed digits.
"""

from __future__ import division, print_function
import numpy as np
from ._constants import eps
from ._utils import empty, mfloor


def bquantize(x, nsd=3, abstol=eps, reltol=10 * eps):
    """Bidirectionally quantize a 1D vector ``x`` to ``nsd`` signed digits.

    This method will terminate early if the error is less than the specified
    tolerances.

    The quantizer details are repeated here for the user's convenience:

        The quantizer is ideal, producing integer outputs centered about zero.
        Quantizers with an even number of levels are of the mid-rise type and
        produce outputs which are odd integers. Quantizers with an odd number
        of levels are of the mid-tread type and produce outputs which are even
        integers.

        .. image:: ../doc/_static/quantizer_model.png
            :align: center
            :alt: Quantizer model

    **Parameters:**

    x : array_like or sequence
           the data to be quantized.

    nsd : int, optional
          The number of signed digits.

    abstol and reltol : floats, optional
        If not supplied, the absolute tolerance and the relative 
        tolerance default to ``eps`` and ``10*eps``, resp.

    **Returns:**

    y : list 
        List of objects described below.

    ``y`` is a list of instances with the same length as ``x`` and the 
    following attributes:

    * ``y[i].val`` is the quantized value in floating-point form,
    * ``y[i].csd`` is a 2-by-nsd (or less) matrix containing
      the powers of two (first row) and their signs (second row).
    
    .. seealso:: 
        :func:`bunquantize`, :func:`ds_quantize`

    """

    n = x.shape[0] if isinstance(x, np.ndarray) else len(x)
    #q = np.zeros((2*n, nsd)) in the original source #rep?
    y = [empty() for i in range(n)]
    offset = -np.log2(0.75)

    for i in range(n):
        xp = x[i]
        y[i].val = 0.
        y[i].csd = np.zeros((2, 0), dtype='int16')
        for _ in range(nsd):
            error = np.abs(y[i].val - x[i])
            if error <= abstol and error <= np.abs(x[i]) * reltol:  # rep? in the orig: or
                break
            p = mfloor(np.log2(np.abs(xp)) + offset)
            p2 = 2 ** p
            sx = np.sign(xp)
            xp = xp - sx * p2
            y[i].val = y[i].val + sx * p2
            addme = np.array((p, sx)).reshape((2, 1))
            y[i].csd = np.concatenate((y[i].csd, addme), axis=1)
    return y


def test_bquantize():
    """Test function for bquantize()
    """
    import scipy.io
    import pkg_resources
    x = np.linspace(-10, 10, 101)
    y = bquantize(x)
    yval = [yi.val for yi in y]
    ycsd = [yi.csd for yi in y]
    fname = pkg_resources.resource_filename(__name__, "test_data/test_bquantize.mat")
    s = scipy.io.loadmat(fname)['s']
    mval = []
    mcsd = []
    for i in range(s.shape[1]):
        mval.append(float(s[0, i][0]))
        mcsd.append(s[0, i][1])
    for i in range(len(mval)):
        assert np.allclose(mval[i], yval[i], atol=1e-8, rtol=1e-5)
        print(mcsd[i].shape, ycsd[i].shape)
        assert np.prod(mcsd[i].shape) + np.prod(ycsd[i].shape) == 0 or \
               mcsd[i].shape == ycsd[i].shape

        if 0 not in ycsd[i].shape:
            assert np.allclose(mcsd[i], ycsd[i], atol=1e-8, rtol=1e-5)
