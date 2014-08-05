# -*- coding: utf-8 -*-
# _axisLabels.py
# Module providing the axisLabel function
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

"""Module providing the axisLabel() function
"""

import collections

import numpy as np


def axisLabels(ran, incr):
    """Utility function to quickly generate the alphanum. axis labels.

    **Parameters:**

    ran : sequence
        Sequence containing the axis points (floats)

    incr : int, or 2-elements sequence
        This parameter may be:

    * an int, the function returns an array of strings corresponding to:
      each element of ``range[0]:range[-1]:incr`` formatted as ``'%g'``.

    * a list, the function returns an array of strings corresponding to:
      each element of ``incr[1]:range[-1]:incr[0]`` formatted as ``'%g'``.

    .. note:: All elements in ``ran`` less than 1e-6 are rounded down to 0.

    **Returns:**

    labels : list of strings

    **Raises:**

    ValueError: "Unrecognised incr."
    """
    ran = np.asarray(ran)
    ran[np.abs(ran) < 1e-6] = 0
    s = []
    if not isinstance(incr, collections.Iterable):
        incr = int(incr)
        first = 0
    elif len(incr) == 2:
        first = incr[1]
        incr = incr[0]
    else:
        raise ValueError("Unrecognised incr: " + str(incr))
    for i in range(first, len(ran), incr):
        s += ['%g' % ran[i]]
    return s

