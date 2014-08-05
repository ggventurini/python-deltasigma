# -*- coding: utf-8 -*-
# _l1norm.py
# Module providing the l1norm function
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

"""Module providing the l1norm() function
"""

from __future__ import division

import numpy as np

from warnings import warn

from scipy.signal import lti, dimpulse

from ._utils import _is_zpk, _is_A_B_C_D, _is_num_den, _get_zpk

def l1norm(H):
    """Compute the l1-norm of a z-domain transfer function.

        The norm is evaluated over the first 100 samples.

        **Parameters:**

        H : sequence or lti object
            Any supported LTI representation is accepted.

        **Returns:**

        l1normH : float
            The L1 norm of ``H``.

        .. note:
            LTI objects are translated to ZPK tuples, with possible
            rounding errors.

    """
    if _is_zpk(H):
        z, p, k = H
        HP = (z, p, k, 1.)
    elif _is_num_den(H):
        num, den = H
        HP = (num, den, 1.)
    elif _is_A_B_C_D(H):
        A, B, C, D = H
        HP = (A, B, C, D, 1.)
    elif isinstance(H, lti):
        warn('l1norm() got an LTI object, translated to zpk form, rounding errors possible.')
        z, p, k = _get_zpk(H)
        HP = (z, p, k, 1.)
    _, y = dimpulse(HP, t=np.arange(100))
    return np.sum(np.abs(y[0]))

