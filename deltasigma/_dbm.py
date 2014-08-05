# -*- coding: utf-8 -*-
# _dbm.py
# This module provides the dbm function.
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

"""This module provides the dbm() function, used to convert an RMS voltage to dBm.
"""

from __future__ import division
import numpy as np

from ._utils import carray, save_input_form, restore_input_form


def dbm(v, R=50):
    """Calculate the dBm equivalent of an RMS voltage ``v``.

    .. math::

        P_{dBm} = 10 \\mathrm{log}_{10}(1000 \\frac{v^2}{R} )

    **Parameters:**

    v : scalar or sequence
        The voltages to be converted.

    R : scalar, optional
        The resistor value the power is calculated upon, defaults to 50 ohm.

    **Returns:**

    PdBm : scalar or sequence
           The input in dBm.

    .. seealso:: :func:`undbm`, :func:`db`, :func:`dbp`, :func:`dbv`

    """
    iform = save_input_form(v)
    v = carray(v)
    y = -np.Inf * np.ones(np.size(v))
    nonzero = (v != 0)
    y[nonzero] = 10. * np.log10(np.abs(v[nonzero] ** 2.) / R) + 30
    return restore_input_form(y, iform)

