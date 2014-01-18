# -*- coding: utf-8 -*-
# _thermometer.py
# This module provides the thermometer function.
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

"""This module provides the thermometer() function.
"""

import numpy as np

def thermometer(x, m):
	""" t = thermometer(x,m)	
	t is an m by length(x) matrix wherein the first
	x(i) components of column i are one
	"""
	t = np.zeros((m, len(x)))
	for i in range(len(x)):
	    t[:x[i], i] = np.ones((x[i], ))
	return t
	
def test():
	"""Test function for thermometer()"""
	tv = np.arange(50)
	rm = np.zeros((70, 50))
	for i in range(50):
		rm[:i, i] = np.ones(rm[:i, i].shape)
	assert np.allclose(thermometer(tv, 70), rm, rtol=1e-05, atol=1e-08)

