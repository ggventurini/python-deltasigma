# -*- coding: utf-8 -*-
# _evalF1.py
# Module providing the evalF1 function
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

"""Module providing the evalF1() function
"""
from __future__ import division

import numpy as np

def evalF1(f1, z, phi=None):
    """Calculate the values of the F1 filter 
    (tranformed prototype filter) of a Saramaki HBF at the given points.
    """
    if phi is not None:
        z = z/phi
    f1 = np.asarray(f1).squeeze()
    f1 = np.atleast_1d(f1)

    F1 = 0.5
    for i in range(f1.shape[0]):
        F1 = F1 + f1[i]*z**(2*i+1)

    return F1
