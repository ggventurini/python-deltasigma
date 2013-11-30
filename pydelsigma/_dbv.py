# -*- coding: utf-8 -*-
# _dbv.py
# This module provides the dbv function.
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

"""This module provides the dbv() function, used to convert a voltage gain to dB.
"""

import numpy as np
from ._undbv import undbv


def dbv(x):
	""" dbv(x) = 20*log10(abs(x)); the dB equivalent of the voltage x"""
	if not hasattr(x, 'shape'):
		if not hasattr(x, '__len__'):
			x = np.array((x,))
		else:
			x = np.array(x)
	elif x.shape == ():
		x = np.array((x,))
	y = -np.inf*np.ones(x.shape)
	nonzero = (x != 0)
	y[nonzero] = 20.*np.log10(np.abs(x[nonzero]))
	return y

def test_dbv():
	"""Test function for dbv()"""
	t = np.array([3.0])
	r1 = undbv(dbv(t))
	assert np.allclose(t, r1, atol=1e-8, rtol=1e-5)

