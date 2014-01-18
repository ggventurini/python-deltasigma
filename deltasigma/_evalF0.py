# -*- coding: utf-8 -*-
# _evalF0.py
# Module providing the evalF0 function
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

"""Module providing the evalF0() function
"""

from __future__ import division
from ._evalF1 import evalF1

def evalF0(f1, z, phi):
    """Calculate the values of the F0 (prototype) filter 
    of a Saramaki HBF at the given points.
    """
    return evalF1(f1, 0.5*(z + 1./z), phi)
