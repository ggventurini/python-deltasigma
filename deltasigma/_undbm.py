# -*- coding: utf-8 -*-
# _undbm.py
# This module provides the undbm function.
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

"""This module provides the undbm() function.
"""

import numpy as np

from ._utils import carray, save_input_form, restore_input_form

def undbm(p, z=50):
    """Calculate the RMS voltage equivalent of a power ``p`` expressed in dBm.

    .. math::

        V_{\\mathrm{RMS}} = \\sqrt{z\\ 10^{p/10 - 3}}

    **Parameters:**

    p : scalar or sequence
        The power to be converted.

    z : scalar, optional
        The normalization resistance value, defaults to 50 ohm.

    **Returns:**

    Vrms : scalar or sequence
           The RMS voltage corresponding to p

    .. seealso:: :func:`undbp`, :func:`undbv`, :func:`dbm`, :func:`db`

    """
    iform = save_input_form(p)
    p = carray(p)
    up = np.sqrt(z*10.**(p/10.-3))
    return restore_input_form(up, iform)

