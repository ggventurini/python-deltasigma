# -*- coding: utf-8 -*-
# _dbp.py
# This module provides the dbp function.
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

"""This module provides the dbp() function, used to convert a power gain to dB.
"""

import numpy as np

def dbp(x):
	""" dbp(x) = 10*log10(x): the dB equivalent of the power x"""
	if not hasattr(x, 'shape'):
		if not hasattr(x, '__len__'):
			x = np.array((x,))
		else:
			x = np.array(x)
	elif x.shape == ():
		x = np.array((x,))
	y = -np.inf*np.ones(x.shape)
	nonzero = (x != 0)
	y[nonzero] = 10.*np.log10(np.abs(x[nonzero]))
	return y

def test_dbp():
	"""Test function for dbp()
	"""
	tv = np.array([2])
	r = np.array([3.01029996])
	res = dbp(tv)
	assert np.allclose(r, res, atol=1e-8, rtol=1e-5)
	tv = 2
	r = 3.01029996
	res = dbp(tv)
	assert np.allclose(r, res, atol=1e-8, rtol=1e-5)
	tv = 2, 2
	r = 3.01029996, 3.01029996
	res = dbp(tv)
	assert np.allclose(r, res, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_dbp()
