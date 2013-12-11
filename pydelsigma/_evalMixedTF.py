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

def evalMixedTF(tf, f, df=1e-5):
    """Compute the mixed transfer function tf at a frequency f.
    tf is a struct array with fields Hs and Hz wich represent 
    a continuous-time/discrete-time tfs which must be multiplied together
    and then added up.
    """
    H = np.zeros(f.shape)
    for i in range(max(tf.shape)):
        H = H + evalTFP(tf[i].Hs, tf[i].Hz, f)
    err = np.logical_or(np.isnan(H), np.isinf(H))
    if np.any(err):
        for i in np.where(err):
            H[i] = evalMixedTF(tf, f[i] + df, df*10)
    return H
