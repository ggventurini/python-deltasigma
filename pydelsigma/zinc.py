# -*- coding: utf-8 -*-
# zinc.py
# The zinc function.
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

"""This module provides the zinc() function which calculates the magnitude 
response of a cascade of comb filters.
"""

import numpy as np

def zinc(f, m=64, n=1):
	""" mag = zinc(f, m=64, n=1)	
	Calculate the magnitude response of a cascade of n mth-order 
	comb filters at frequencies f.
	"""
	return np.fabs(np.sinc(m*f)/np.sinc(f) )**n
