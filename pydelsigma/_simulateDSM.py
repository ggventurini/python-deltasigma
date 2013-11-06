# -*- coding: utf-8 -*-
# _simulateDSM.py
# Module providing the simulateDSM function
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

"""Module providing the simulateDSM() function
"""

import numpy as np
from scipy.signal import tf2zpk, zpk2ss
from scipy.linalg import orth, norm, inv
from ._utils import carray

def simulateDSM(u, arg2, nlev=2, x0=0):
	"""[v, xn, xmax, y] = simulateDSM(u, ABCD, nlev=2, x0=0)
	or
	[v, xn, xmax, y] = simulateDSM(u, ntf, nlev=2, x0=0)

	Compute the output of a general delta-sigma modulator with input u,
	a structure described by ABCD, an initial state x0 (default zero) and
	a quantizer with a number of levels specified by nlev.
	Multiple quantizers are implied by making nlev an array,
	and multiple inputs are implied by the number of rows in u.

	Alternatively, the modulator may be described by an NTF.
	The NTF is zpk object. (The STF is assumed to be 1.)
	The structure that is simulated is the block-diagional structure used by
	zp2ss.m.
	"""

	#fprintf(1,'Warning: You are running the non-mex version of simulateDSM.\n');
	#fprintf(1,'Please compile the mex version with "mex simulateDSM.c"\n');

	nlev = carray(nlev)
	u = np.array(u) if not hasattr(u, 'ndim') else u
	if not max(u.shape) == np.prod(u.shape):
		raise ValueErrror("The u vector has shape %s" % u.shape)
	if u.ndim == 1:
		u = u.reshape((1, -1))
	nu = u.shape[0]
	nq = 1 if np.isscalar(nlev) else nlev.shape[0]
	# extract poles and zeros
	if (hasattr(arg2, 'inputs') and not arg2.inputs == 1) or \
	   (hasattr(arg2, 'outputs') and not arg2.outputs == 1):
			raise TypeError("The supplied TF isn't a SISO transfer function.")
	if hasattr(arg2, 'num') and hasattr(arg2, 'den'):
		filt = hasattr(arg2, 'outputs')
		num = arg2.num[0][0] if filt else arg2.num
		den = arg2.den[0][0] if filt else arg2.den
		zeros, poles, k = tf2zpk(num, den)
	elif (hasattr(arg2, 'zeros') and hasattr(arg2, 'poles')) or \
	   (hasattr(arg2, 'zero') and hasattr(arg2, 'pole')):
		# LTI objects have poles and zeros, 
		# TransferFunction-s have pole() and zero()
	   	zeros = arg2.zeros if hasattr(arg2, 'zeros') else arg2.zero()
	   	poles = arg2.poles if hasattr(arg2, 'poles') else arg2.pole()
		if hasattr(arg2, 'k'):
			k = arg2.k
		elif hasattr(arg2, 'gain'):
			k = arg2.gain  
		elif hasattr(arg2, 'returnScipySignalLti'): 
			k = np.array(arg2.returnScipySignalLti()[0][0].gain)
	elif hasattr(arg2, 'form') and arg2.form == 'coeff':
		zeros, poles, k = tf2zpk(arg2.num, arg2.den)
	elif hasattr(arg2, 'form'):
		raise ValueError('%s: Unknown form: %s' % (__name__, arg2.form))
	elif hasattr(arg2, '__len__'):
		if len(arg2) == 3:
			zeros, poles, k = arg2
		else:
			ABCD = carray(arg2)
			if ABCD.shape[1] != ABCD.shape[0] + nu:
				raise ValueError('The ABCD argument does not have proper dimensions.')
	else:
		raise TypeError('%s: Unknown transfer function %s' % (__name__, str(arg2)))
		
	# need to set order and form now.
	form = 2 - 1*(not 'zeros' in locals())
	order = carray(zeros).shape[0] if form == 2 else ABCD.shape[0] - nq
	
	if not hasattr(x0 , '__len__'):
		x0 = x0*np.ones((order, 1))
	else:
		x0 = np.array(x0).reshape((-1, 1))
	
	if form == 1:
		A = ABCD[:order, :order]
		B = ABCD[:order, order:order+nu+nq]
		C = ABCD[order:order+nq, :order]
		D1 = ABCD[order:order+nq, order:order+nu]
	else:
		A, B2, C, D2 = zpk2ss(poles, zeros, -1)	# A realization of 1/H
		# Transform the realization so that C = [1 0 0 ...]
		Sinv = orth(np.hstack((np.transpose(C), np.eye(order)))) / norm(C)
		S = inv(Sinv)
		C = np.dot(C, Sinv)
		if C[0, 0] < 0:
			S = -S
			Sinv = -Sinv
		A = np.dot(np.dot(S, A), Sinv) 
		B2 = np.dot(S, B2) 
		C = np.hstack((np.ones((1, 1)), np.zeros((1, order-1)))) # C=C*Sinv; 
		D2 = np.zeros((0,))
		# !!!! Assume stf=1
		B1 = -B2
		D1 = 1
		B = np.hstack((B1, B2))

	N = u.shape[1]
	v = np.empty((nq, N))
	y = np.empty((nq, N)) 	# to store the quantizer input
	xn = np.empty((order, N), dtype=np.complex64) # to store the state information 
	xmax = np.abs(x0) # to keep track of the state maxima

	for i in range(N):
		# y0 needs to be cast to real because ds_quantize needs real
		# inputs. If quantization were defined for complex numbers,
		# this cast could be removed
		y0 = np.real(np.dot(C, x0) + np.dot(D1, u[:, i]))
		y[:, i] = y0
		v[:, i] = ds_quantize(y0, nlev)
		x0 = np.dot(A, x0) + np.dot(B, np.vstack((u[:, i], v[:, i])))
		xn[:, i] = x0.T
		xmax = np.max((np.abs(x0), xmax), 0)

	return v.squeeze(), xn.squeeze(), xmax, y.squeeze()

def ds_quantize(y, n):
	"""v = ds_quantize(y,n)
	Quantize y to:
	 
	* an odd integer in [-n+1, n-1], if n is even, or
	* an even integer in [-n, n], if n is odd.

	This definition gives the same step height for both mid-rise
	and mid-tread quantizers.
	"""
	v = np.zeros(y.shape)
	for qi in range(n.shape[0]): 
		if n[qi] % 2 == 0: # mid-rise quantizer
			v[qi, 0] = 2*np.floor(0.5*y[qi, 0]) + 1
		else: # mid-tread quantizer
			v[qi, 0] = 2*np.floor(0.5*(y[qi, 0] + 1))
		L = n[qi] - 1
		v[qi, 0] = np.sign(v[qi, 0])*np.max((np.abs(v[qi, 0]), L))
	return v
