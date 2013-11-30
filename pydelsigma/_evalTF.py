# -*- coding: utf-8 -*-
# _evalTF.py
# This module provides the evalTF function.
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

"""This module provides the evalTF() function, used to evaluate a tf at the 
point(s) given by the user.
"""

import numpy as np
from scipy.signal import tf2zpk
from ._evalRPoly import evalRPoly

def evalTF(tf, z):
	"""h = evalTF(tf, z)
	Evaluates the rational function described by tf at the point(s) given in the z vector.
	
	TF must be either a TransferFunction or an LTI object or any object having the attributes:
	 zeros, poles, k	if form == 'zp'
	or:
	 num, den		if form == 'coeff'
	"""
	# Original comment in deltasig:
	# In Matlab 5, the ss/freqresp() function does nearly the same thing.
	
	# in our case a transfer function is a 'TransferFunction' not a zpk
	if (hasattr(tf, 'inputs') and not tf.inputs == 1) or \
	   (hasattr(tf, 'outputs') and not tf.outputs == 1):
		raise TypeError, "Only SISO transfer functions can be evaluated."
	if hasattr(tf, 'num') and hasattr(tf, 'den'):
		# for now we support both TransferFunction objects (python-control)
		# and lti objects (scipy).
		filt = hasattr(tf, '__class__') and tf.__class__.__name__ == 'TransferFunction'
		num = tf.num[0][0] if filt else tf.num
		den = tf.den[0][0] if filt else tf.den
		h = np.polyval(num, z) / np.polyval(den, z)
	elif (hasattr(tf, 'zeros') and hasattr(tf, 'poles')) or \
	   (hasattr(tf, 'zero') and hasattr(tf, 'pole')):
		# LTI objects have poles and zeros, 
		# TransferFunction-s have pole() and zero()
		zeros = tf.zeros if hasattr(tf, 'zeros') else tf.zero()
		poles = tf.poles if hasattr(tf, 'poles') else tf.pole()
		if hasattr(tf, 'k'):
			k = tf.k
		elif hasattr(tf, 'gain'):
			k = tf.gain  
		elif hasattr(tf, 'returnScipySignalLti'): 
			k = np.array(tf.returnScipySignalLti()[0][0].gain)
		h = k * evalRPoly(zeros, z) / evalRPoly(poles, z)
	elif hasattr(tf, 'form') and tf.form == 'zp':
		h = tf.k * evalRPoly(tf.zeros, z) / evalRPoly(tf.poles, z)
	elif hasattr(tf, 'form') and tf.form == 'coeff':
		h = np.polyval(tf.num, z) / np.polyval(tf.den, z)
	elif hasattr(tf, 'form'):
		raise ValueError, '%s: Unknown form: %s' % (__name__, tf.form)
	elif hasattr(tf, '__len__'):
		if len(tf) == 2:
			num, den = tf[0], tf[1]
			h = np.polyval(num, z) / np.polyval(den, z)
		elif len(tf) == 3:
			zeros, poles, k = tf
			h = k * evalRPoly(zeros, z) / evalRPoly(poles, z)
	else:
		raise TypeError, '%s: Unknown transfer function %s' % (__name__, str(tf))
	return h
	
def test_evalTF():
	"""Test function for evalTF()
	"""
	from scipy.signal import tf2zpk
	from ._utils import empty
	num, den = np.poly([3, 0.3, 1]), np.poly([2, 0.5, .25])
	H = (num, den)
	tstr1 = empty()
	tstr1.form, tstr1.num, tstr1.den = 'coeff', num, den
	tstr2 = empty()
	tstr2.form = 'zp'
	tstr2.zeros, tstr2.poles, tstr2.k = tf2zpk(num, den)
	z = np.exp(1j*np.linspace(0, 2*np.pi, num=129, endpoint=True))
	h1 = evalTF(tstr1, z)
	h2 = evalTF(tstr2, z)
	h3 = evalTF(H, z)
	assert np.allclose(np.abs(h1), np.abs(h2), atol=1e-8, rtol=1e-5)
	assert np.allclose(h1, h3, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_evalTF()
	
