# -*- coding: utf-8 -*-
# _dbp.py
# This module provides the dbp function.
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

"""This module provides the dbp() function, used to convert a power gain to dB.
"""

import numpy as np

from ._utils import carray, save_input_form, restore_input_form


def dbp(x):
    """Calculate the dB equivalent of the power ratio ``x``.

    .. math::

        P_{dB} = 10 \\mathrm{log}_{10}(x)

    **Parameters:**

    x : scalar or sequence
        The power ratio to be converted.

    **Returns:**

    PdB : scalar or sequence
           The input expressed in dB.

    .. seealso:: :func:`undbp`, :func:`db`, :func:`dbm`, :func:`dbv`

    """
    iform = save_input_form(x)
    x = carray(x)
    y = -np.Inf * np.ones(x.shape)
    nonzero = (x != 0)
    y[nonzero] = 10. * np.log10(np.abs(x[nonzero]))
    return restore_input_form(y, iform)

