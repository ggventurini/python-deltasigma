# -*- coding: utf-8 -*-
# _evalMixedTF.py
# Module providing the evalMixedTF function
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

"""Module providing the evalMixedTF() function
"""

from __future__ import division
import numpy as np

from ._evalTFP import evalTFP
from ._utils import carray, save_input_form, restore_input_form

def evalMixedTF(tf, f, df=1e-5):
    """Compute a mixed transfer function.

    Mathematically, it means to evaluate the sum of products:

    .. math::

        TF(f) = \\sum_i H_{z,i}(f) \\cdot H_{s,i}(f)

    **Parameters:**

    tf : dict
        ``tf`` is a dictionary of lists of 1d arrays, with keys 'Hs' and 'Hz',
        which represent continuous-time and discrete-time TFs which will be
        evaluated, multiplied together and then added up.
    f : scalar or sequence-like
        The frequencies (or frequency, if ``f`` is a scalar) at which the
        product will be evaluated.
    df : float, optional
        If the method happens to evaluate a transfer function in a root, it
        will move away of ``df``, defaulting to 1E-5.

    **Returns:**

    TF : scalar or sequence
        The sum of products computed at ``f``.

    """
    iform = save_input_form(f)
    f = carray(f)
    H = np.zeros(f.shape)
    for i in range(max(len(tf['Hs']), len(tf['Hz']))):
        H = H + evalTFP(tf['Hs'][i], tf['Hz'][i], f)
    err = np.logical_or(np.isnan(H), np.isinf(H))
    if np.any(err):
        for i in np.where(err):
            H[i] = evalMixedTF(tf, f[i] + df, df*10)
    return restore_input_form(H, iform)
