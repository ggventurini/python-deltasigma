# -*- coding: utf-8 -*-
# undbm.py
# This module provides the undbm function.
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

"""This module provides the undbm() function.
"""

import numpy as np

def undbm(p, z=50):
	""" v = undbm(p, z=50) = sqrt(z*10^(p/10-3))
	RMS voltage equivalent to a power p in dBm
	z is the normalization resistance (defaults to 50ohm)
	"""
	return np.sqrt(z*10.**(p/10.-3))

def test():
	assert np.allclose([undbm(53.015)], [100.054125892], rtol=1e-05, atol=1e-08)
	assert np.allclose([undbm(3, 100)], [0.44668359215], rtol=1e-05, atol=1e-08)

if __name__ == '__main__':
	test()
