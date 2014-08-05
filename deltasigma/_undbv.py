# -*- coding: utf-8 -*-
# _undbv.py
# This module provides the undbv function.
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

"""This module provides the undbv() function.
"""
from __future__ import division
import numpy as np

from ._utils import carray, save_input_form, restore_input_form

def undbv(x):
    """Convert ``x`` from dB to voltage.

    The conversion is carried out according to the relationship:

    .. math::

        V_{\\mathrm{RMS}} = 10^{x/20}

    **Parameters:**

    x : scalar or sequence
        The signal in dB to be converted.

    **Returns:**

    Vrms : scalar or sequence
           The RMS voltage corresponding to x.

    .. seealso:: :func:`undbm`, :func:`undbp`, :func:`dbv`, :func:`db`

    """
    iform = save_input_form(x)
    x = carray(x)
    up = 10.**(x/20.)
    return restore_input_form(up, iform)
    
