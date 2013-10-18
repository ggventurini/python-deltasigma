# -*- coding: utf-8 -*-
# _utils.py
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

from warnings import warn
import fractions
from fractions import Fraction as Fr
import numpy as np
from ._dbp import dbp
from ._dbv import dbv
from._constants import eps

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
	"""MATLAB-like floor function.
	"""
	def _mfloor(z):
		return np.floor(z) if np.sign(z) >= 0 else -np.ceil(-z)
	_internal = np.frompyfunc(_mfloor, 1, 1)
	return np.array(_internal(x), dtype=x.dtype)

def zpk(z, p, k):
	"""Returns a zpk object with the zeros (z), poles (p) and gain (k) provided.
	"""
	t = empty()
	t.form = 'zp'
	t.zeros, t.poles, t.k = carray(z), carray(p), k
	return t
	
def db(x, input_type='voltage', R=1.):
	"""The return value is defined as:
	y = dbv(x) - 20*log10(R) if input_type == 'voltage'
	y = dbp(x) if input type
	
	MATLAB provides a function with this exact signature.
    """
	if input_type.lower().strip() == 'voltage':
			y = dbv(x) - 10.*np.log10(R)
	elif input_type.lower().strip() == 'power':
		y = dbp(x)
		if R != 1.:
			warn("db called with a non default R value, " + 
			     "but R is going to be ignored since input_type is power",
				 RuntimeWarning)
	else:
		raise ValueError, "db got input_type %s, instead of voltage or power" % input_type
	return y

def carray(x):
	"""Check that x is an ndarray. If not, try to convert it to ndarray.
	"""
	if not hasattr(x, 'shape'):
		if not hasattr(x, '__len__'):
			x = np.array((x,))
		else:
			x = np.array(x)
	elif x.shape == ():
		x = np.array((x,))
	else:
		pass #nothing to do here
	return x

def cplxpair(x, tol=100):
	"""
	Sort complex numbers into complex conjugate pairs.

	This function replaces MATLAB's cplxpair for vectors.
	"""
	x = carray(x)
	x = x.tolist()
	x = map(lambda x: np.real_if_close(x, tol), x)
	xreal = np.array(filter(np.isreal, x))
	xcomplex = np.array(filter(np.iscomplex, x))
	xreal = np.sort_complex(xreal)
	xcomplex = np.sort_complex(xcomplex)
	xcomplex_ipos = xcomplex[xcomplex.imag >  0.]
	xcomplex_ineg = xcomplex[xcomplex.imag <= 0.]
	res = []
	for i, j in zip(xcomplex_ipos, xcomplex_ineg):
		if not abs(i - np.conj(j)) < tol*eps:
			raise ValueError, "Complex numbers can't be paired."
		res += [j, i]
	return np.hstack((np.array(res), xreal))

	
def test_rat():
	"""Test function for rat()"""
	import numpy.random as rnd
	for i in range(10):
		n, d = rnd.randint(1, 5000), rnd.randint(1, 5000)
		fr = float(n)/float(d)
		nt, dt = rat(fr, tol=200e-6)
		assert np.allclose((n/float(d) -  nt/float(dt),), (0.,), atol=200e-6, rtol=1e-12)

def test_gcd_lcm():
	"""Test function for gcd() and lcm"""
	a, b = 36, 31721
	tlcm, tgcd = 1141956, 1
	assert lcm(a, b) == tlcm
	assert gcd(a, b) == tgcd

def test_mfloor():
	"""Test function for mfloor()"""
	tv = np.linspace(-1, 1, 10)
	tres = np.zeros(tv.shape)
	tres[tv>=0] = np.floor(tv[tv>=0])
	tres[tv<0] = -np.ceil(np.abs(tv[tv<0]))
	tresf = mfloor(tv)
	assert np.allclose(tres, tresf, atol=1e-8, rtol=1e-5)

def test_zpk():
	"""Test function for zpk.
	"""
	z = [2,]
	p = [1, 3]
	k = 1
	t = zpk(z, p, k)
	assert t.form == 'zp'
	assert t.zeros.tolist() == z
	assert t.poles.tolist() == p
	assert t.k == k	

def test_db():
	import warnings
	from ._undbv import undbv
	tv = np.array([2])
	r = np.array([3.01029996])
	res = db(tv, 'power')
	assert np.allclose(r, res, atol=1e-8, rtol=1e-5)
	tv = 2
	r = 3.01029996
	a = False
	res = db(tv, 'power') 
	assert np.allclose(r, res, atol=1e-8, rtol=1e-5)
	tv = 2, 2
	r = 3.01029996, 3.01029996
	res = db(tv, 'power')
	assert np.allclose(r, res, atol=1e-8, rtol=1e-5)
	t = np.array([3.0])
	r1 = undbv(db(t, 'voltage'))
	assert np.allclose(t, r1, atol=1e-8, rtol=1e-5)
	
def test_cplxpair():
	a = np.array([1 + eps*20j, 1.1 + 2j, 1.1 - (2+50*eps)*1j, .1 + (1+99*eps)*.2j, .1 - .2j])
	assert np.allclose(cplxpair(a), np.array([0.1-0.2j, 0.1+0.2j, 1.1-2.j, 1.1+2.j, 1.0+0.j]), 
	                   atol=100*eps)

if __name__ == '__main__':
	test_rat()
	test_gcd_lcm()
	test_mfloor()
	test_zpk()
	test_db()
	test_cplxpair()
