# -*- coding: utf-8 -*-
# _evalTF.py
# This module provides the evalTF function.
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

"""This module provides the evalTF() function, used to evaluate a tf at the 
point(s) given by the user.
"""

import numpy as np
from ._evalRPoly import evalRPoly
from ._utils import _is_num_den, _is_zpk, _get_zpk

def evalTF(tf, z):
    """Evaluates the rational function ``tf`` at the point(s)
    given in ``z``.

    **Parameters:**

    tf : object
        the LTI description of the CT system, which can be in one of the
        following forms:

    * an LTI object,
    * a zpk tuple,
    * a (num, den) tuple,
    * an ABCD matrix (internally converted to zpk representation),
    * a list-like containing the A, B, C, D matrices (also internally converted to zpk representation).

    z : scalar or 1d ndarray
        The z values for which ``tf`` is to be evaluated.

    **Returns:**

    tf(z) : scalar or 1d-ndarray
        The result.

    """
    # Original comment in deltasig:
    # In Matlab 5, the ss/freqresp() function does nearly the same thing.

    # in our case a transfer function is a scipy LTI object
    if (hasattr(tf, 'inputs') and not tf.inputs == 1) or \
       (hasattr(tf, 'outputs') and not tf.outputs == 1):
        raise TypeError("Only SISO transfer functions can be evaluated.")

    if hasattr(tf, 'num') and hasattr(tf, 'den'):
        h = np.polyval(tf.num, z) / np.polyval(tf.den, z)
    elif (hasattr(tf, 'zeros') and hasattr(tf, 'poles')):
        # LTI objects have poles and zeros
        zeros = tf.zeros
        poles = tf.poles
        k = tf.gain
        h = k*evalRPoly(zeros, z)/evalRPoly(poles, z)
    # we try not to convert the given representation to another one
    elif _is_num_den(tf):
        num, den = tf[0], tf[1]
        h = np.polyval(num, z)/np.polyval(den, z)
    elif _is_zpk(tf):
        zeros, poles, k = tf
        h = k*evalRPoly(zeros, z)/evalRPoly(poles, z)
    else:
        # ABCD and A, B, C, D will be converted through _utils
        zeros, poles, k = _get_zpk(tf)
        h = k*evalRPoly(zeros, z)/evalRPoly(poles, z)
    return h

