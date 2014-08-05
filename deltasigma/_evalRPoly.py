# -*- coding: utf-8 -*-
# _evalRPoly.py
# This module provides the evalRPoly function.
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

"""This module provides the evalRPoly() function, used to evaluate the value 
of a polynomial which is given in terms of its roots.
"""

import numpy as np
from ._utils import carray

def evalRPoly(roots, x, k=1):
	"""Compute the value of a polynomial which is given in terms of its roots.
	"""
	roots = carray(roots)
	y = k
	roots = roots[~np.isinf(roots)]        # remove roots at infinity
	for r in roots:
		y = y*(x - r)
	return y

