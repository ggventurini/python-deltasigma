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
from scipy.signal import lti
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
	"""An empty function used to hold attributes"""
	def __init__(self):
		pass

def mfloor(x):
	"""MATLAB-like floor function.
	"""
	def _mfloor(z):
		"""Base function to generate the ufunc floor"""
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

def minreal(tf, tol=None):
	"""Remove pole/zero pairs from a transfer function
	when the two match within the tolerance tol.
	
	tf may be a transfer function (lti object) or a list of transfer functions
	tol is optional and defaults to the system epsilon if unset.
	
	return a list of tfs or a tf, depending on the input type.
	"""
	# initially based on python-control
	# which is in turn based on octave minreal
	# then modified considerably

	# recursively handle multiple tfs
	if isinstance(tf, list) or isinstance(tf, tuple):
		ret = []
		for tfi in tf:
			ret += [minreal(tfi, tol)]
		return ret

	# default accuracy
	sqrt_eps = np.sqrt(eps)
	
	if (hasattr(tf, 'inputs') and not tf.inputs == 1) or \
	   (hasattr(tf, 'outputs') and not tf.outputs == 1):
		raise TypeError, "Only SISO transfer functions can be evaluated."
	if hasattr(tf, 'num') and hasattr(tf, 'den'):
		#filt = hasattr(tf, 'outputs')
		num = tf.num #[0][0] if filt else tf.num
		den = tf.den #[0][0] if filt else tf.den
		k = num[0] / den[0]
		zeros = np.roots(num)
		poles = np.roots(den)
	elif (hasattr(tf, 'zeros') and hasattr(tf, 'poles')) or \
	   (hasattr(tf, 'zero') and hasattr(tf, 'pole')):
		# LTI objects have poles and zeros, 
		# python-control TransferFunction-s have pole() and zero()
		zeros = tf.zeros if hasattr(tf, 'zeros') else tf.zero()
		poles = tf.poles if hasattr(tf, 'poles') else tf.pole()
		if hasattr(tf, 'k'):
			k = tf.k
		elif hasattr(tf, 'gain'):
			k = tf.gain  
		elif hasattr(tf, 'returnScipySignalLti'): 
			k = np.array(tf.returnScipySignalLti()[0][0].gain)
	else:
		raise ValueError, "Unknown transfer function type."
	
	zeros = carray(zeros)
	poles = carray(poles)
	zeros.sort()
	poles.sort()
	reducedzeros = []

	for z in zeros:
		t = tol or 1000 * max(eps, abs(z) * sqrt_eps)
		idx = np.where(abs(poles - z) < t)[0]
		if len(idx):
			# cancel this zero against one of the poles
			# remove the pole and do not add the zero to the new
			poles = np.delete(poles, idx[0])
		else:
			# no matching pole
			reducedzeros.append(z)
	if len(reducedzeros):
		newzeros = carray(reducedzeros)
		num = k * np.real(np.poly(newzeros))
	else:
		num = np.array([k])
	den = np.real(np.poly(poles))

	return lti(num, den)

def diagonal_indices(a, offset=0):
	"""Return the indices to the main diagonal of a 2D array a (if offset = 0),
	or to a secondary diagonal, having the offset from the main one as specified.

	The array a does not need to be square.

	Note: the sup-diagonal is at offset +1, the sub-diagonal at offset -1.

	Returns: (xs, ys)
	"""
	di, dj = np.diag_indices_from(a[:min(a.shape), :min(a.shape)])
	if offset > 0:
		di, dj = zip(*[(i, j) for i, j in zip(di, dj+offset) if 0 <= j < a.shape[1]])
	elif offset == 0:
		pass
	else:
		di, dj = zip(*[(i, j) for i, j in zip(di-offset, dj) if 0 <= i < a.shape[0]])
	return di, dj

def circshift(a, shift):
	"""Shift array circularly.
	
    The circshift(a, shift) function circularly shifts the values in the 
	array 'a' by 'shift' elements. 
	
	Parameters:
	
	a, ndarray
	       the array to be shifted. Notice that a should have a greater or equal 
		   number of dimensions than 'shift' ('shift' being a scalar is equal to 
		   'shift' having one dimension.)
	
	shift, int or ndarray-like of int.
	       the N-th element specifies the shift amount for the N-th dimension 
		   of the input array 'a'. 
		   If an element in 'shift' is positive, the values of A are
           shifted to higher-index rows (ie down) or to higher-index columns 
		   (ie to the right). 
		   If the element is negative, the values of A are shifted in the opposite 
		   directions, towards lower-index rows (ie up) or to lower-index columns 
		   (ie right).
		   If shift is an int, the shift happens along axis=0.
		   All dimensions that do not have a corresponding shift value in 'shift' 
		   are left untouched (ie shift=(1,0,0) is equal to shift=(1,), with the
		   exception that the former will trigger an IndexError if a.ndim < 3)

	Returns:

	The shifted array.
	"""
	if np.isscalar(shift):
		shift = [shift]
	for axis, ashift in enumerate(shift):
		a = np.roll(a, ashift, axis=axis)
	return a

def test_rat():
	"""Test function for rat()"""
	import numpy.random as rnd
	for _ in range(10):
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
	tres[tv >= 0] = np.floor(tv[tv>=0])
	tres[tv <  0] = -np.ceil(np.abs(tv[tv<0]))
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
	"""Test function for db()
	"""
	from ._undbv import undbv
	tv = np.array([2])
	r = np.array([3.01029996])
	res = db(tv, 'power')
	assert np.allclose(r, res, atol=1e-8, rtol=1e-5)
	tv = 2
	r = 3.01029996
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
	"""Test function for cplxpair()
	"""
	a = np.array([1 + eps*20j, 1.1 + 2j, 1.1 - (2+50*eps)*1j, .1 + (1+99*eps)*.2j, .1 - .2j])
	assert np.allclose(cplxpair(a), np.array([0.1-0.2j, 0.1+0.2j, 1.1-2.j, 1.1+2.j, 1.0+0.j]), 
	                   atol=100*eps)

def test_diagonal_indices():
	"""Unit test for diagonal_indices()
	"""
	a = np.arange(1, 26)
	a = a.reshape((5, 5))
	d = a[diagonal_indices(a)]
	dp1 = a[diagonal_indices(a, 1)]
	dm1 = a[diagonal_indices(a, -1)]
	assert np.allclose(d, (1, 7, 13, 19, 25))
	assert np.allclose(dp1, (2, 8, 14, 20))
	assert np.allclose(dm1, (6, 12, 18, 24))
	a = np.arange(1, 21)
	a = a.reshape((4, 5))
	d = a[diagonal_indices(a)]
	dp1 = a[diagonal_indices(a, 1)]
	dm1 = a[diagonal_indices(a, -1)]
	assert np.allclose(d, (1, 7, 13, 19))
	assert np.allclose(dp1, (2, 8, 14, 20))
	assert np.allclose(dm1, (6, 12, 18))
	a = a.reshape((5, 4))
	d = a[diagonal_indices(a)]
	dp1 = a[diagonal_indices(a, 1)]
	dm1 = a[diagonal_indices(a, -1)]
	assert np.allclose(d, (1, 6, 11, 16))
	assert np.allclose(dp1, (2, 7, 12))
	assert np.allclose(dm1, (5, 10, 15, 20))
	
def test_circshift():
	A = np.arange(1, 10).reshape((3, 3))
	shift = np.array((1, -1))
	Ashifted = circshift(A, shift)
	Ares = np.array(((8, 9, 7), (2, 3, 1), (5, 6, 4)))
	assert np.allclose(Ashifted, Ares)