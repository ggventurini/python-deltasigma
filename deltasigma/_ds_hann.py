# -*- coding: utf-8 -*-
# _ds_hann.py
# This module provides the ds_hann function.
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

"""This module provides the ds_hann() function, used to generate a Hann 
window that does not smear tones located exactly in a bin.
"""

import numpy as np

def ds_hann(n):
    """A Hann window of length :math:`n`. 

    The Hann window, aka *the raised cosine window*, is defined as:

    .. math::

        w(x) = 0.5\\ \\left(1 - cos\\left(\\frac{2 \\pi x}{n}\\right) \\right)

    This windowing function does not smear tones located exactly in a bin.

    **Parameters:**

    n : integer
        The window length, in number of samples.

    **Returns:**

    w : 1d nd_array
        The Hann window.

    .. note:: Functionally equivalent to numpy's ``hanning()``, provided
        to ease porting of code from MATLAB. Also, we take care always to
        return an array of dimensions ``(n,)`` and type ``float_``.

    .. plot::

      import pylab as plt
      from deltasigma import ds_hann
      x = ds_hann(100)
      plt.figure(figsize=(12, 5))
      plt.plot(x, 'o-')
      ax = plt.gca()
      ax.set_ylim(0, 1.02)
      plt.grid(True)
      plt.title("100-samples Hann window")
      plt.xlabel("Sample #")
      plt.ylabel("Value")

    """
    x = np.arange(n, dtype='float_')
    return .5*(1 - np.cos(2*np.pi*x/n))

