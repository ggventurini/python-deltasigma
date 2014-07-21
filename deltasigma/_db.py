# -*- coding: utf-8 -*-
# _db.py
# This module provides the db function.
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

"""This module provides the db() function, used to convert an RMS voltage,
expressed in Volt, or a power, expressed in Watt, to dBm.
"""

from __future__ import division
from warnings import warn
import numpy as np

from ._dbp import dbp
from ._dbv import dbv


def db(x, input_type='voltage', R=1.):
    """Calculate the dB equivalent of the RMS signal ``x``.

    For input type ``"voltage"``, the return value is defined as

    .. math::

        P_{dB} = 20 \\mathrm{log}_{10}\\left(\\frac{x}{R}\\right)


    Otherwise, for input type ``"power"``,

    .. math::

        P_{dB} = 10 \\mathrm{log}_{10}(x)


    **Parameters:**

    x : scalar or sequence
        The signal to be converted.

    input_type : string, optional
        The input type, either ``"voltage"`` or ``"power"``.

    R : float, optional
        The normalization resistor value, used only for voltage inputs.


    **Returns:**

    PdB : scalar or sequence
        The input expressed in dB.


    .. note:: MATLAB provides a function with this exact signature.

    .. seealso:: :func:`undbm`, :func:`undbv`, :func:`undbp`, :func:`dbv`, :func:`dbp`, :func:`dbv`

    """
    if input_type.lower().strip() == 'voltage':
        y = dbv(x) - 10. * np.log10(R)
    elif input_type.lower().strip() == 'power':
        y = dbp(x)
        if R != 1.:
            warn("db called with a non default R value, " +
                 "but R is going to be ignored since input_type is power",
                 RuntimeWarning)
    else:
        raise ValueError("db got input_type %s, instead of voltage or power" % input_type)
    return y

