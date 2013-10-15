# -*- coding: utf-8 -*-
# utils.py
# Miscellaneous functions and stdlib wrappers for MATLAB functions 
# that do not find a direct replacement in numpy/scipy.
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

""" Miscellaneous functions and wrappers for MATLAB functions 
 that do not find a direct replacement in numpy/scipy.
"""

import numpy as np
import fractions
from fractions import Fraction as Fr

def rat(x, tol):
	"""num, den = rat(x, tol)
	where num/den == x to the specified tolerance tol
	Note: num, den are of type 'int'
	"""
	return Fr(float(x)).limit_denominator(int(1/float(tol))).numerator, \
                Fr(float(x)).limit_denominator(int(1/float(tol))).denominator

gcd = fractions.gcd

lcm = lambda a, b: int(a*b / float(gcd(a, b)))

class empty:
	pass

def mfloor(x):
	def _mfloor(z):
		return np.floor(z) if np.sign(z) >= 0 else -np.ceil(-z)
	_internal = np.frompyfunc(_mfloor, 1, 1)
	return np.array(_internal(x), dtype=x.dtype)


def test_rat():
	import numpy.random as rnd
	for i in range(10):
		n, d = rnd.randint(1, 5000), rnd.randint(1, 5000)
		fr = float(n)/float(d)
		nt, dt = rat(fr, tol=200e-6)
		assert np.allclose((n/float(d) -  nt/float(dt),), (0.,), atol=200e-6, rtol=1e-12)

def test_gcd_lcm():
	a, b = 36, 31721
	tlcm, tgcd = 1141956, 1
	assert lcm(a, b) == tlcm
	assert gcd(a, b) == tgcd

def test_mfloor():
	tv = np.linspace(-1, 1, 10)
	tres = np.zeros(tv.shape)
	tres[tv>=0] = np.floor(tv[tv>=0])
	tres[tv<0] = -np.ceil(np.abs(tv[tv<0]))
	tresf = mfloor(tv)
	assert np.allclose(tres, tresf, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_rat()
	test_gcd_lcm()
	test_mfloor()
