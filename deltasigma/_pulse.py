# -*- coding: utf-8 -*-
# _pulse.py
# This module provides the pulse function.
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

"""This module provides the pulse() function, which calculates the sampled 
pulse response of a CT system.
"""

from __future__ import division
import numpy as np
from scipy.signal import step2, lti

from ._utils import lcm, rat

def pulse(S, tp=(0., 1.), dt=1., tfinal=10., nosum=False):
	""" y = pulse(S, tp=[0 1], dt=1, tfinal=10, nosum=0)
	Calculate the sampled pulse response of a CT system. 
	tp may be an array of pulse timings, one for each input.
	Outputs:
	y	The pulse response
	
	Inputs
	S	A **list** of LTI objects specifying the system.
	tp	An (n, 2) array of pulse timings
	dt	The time increment
	tfinal	The time of the last desired sample
	nosum	A flag indicating that the responses are not to be summed
	"""
	tp = np.mat(np.array(tp))
	if len(tp.shape) == 1:
		if not tp.shape[0] == 2:
			raise ValueError("tp is not (n, 2)-shaped")
		tp.reshape((1, tp.shape[0]))
	if len(tp.shape) == 2: 
		if not tp.shape[1] == 2:
			raise ValueError("tp is not (n, 2)-shaped")

	# Compute the time increment
	dd = 1;
	for tpi in np.nditer(tp.T.copy(order='C')):
		x, di = rat(tpi, 1e-3)
		dd = lcm(di, dd)

	x, ddt = rat(dt, 1e-3)
	x, df = rat(tfinal, 1e-3)
	delta_t = 1. / lcm(dd, lcm(ddt, df))
	delta_t = max(1e-3, delta_t)	# Put a lower limit on delta_t
	y1 = None
	for Si in S:
		T1, y1i= step2(Si, T=np.arange(0., tfinal + delta_t, delta_t))
		if y1 is None:
			y1 = y1i.reshape((y1i.shape[0], 1, 1))
		else:
			y1 = np.concatenate((y1, 
                                             y1i.reshape((y1i.shape[0], 1, 1))), 
                                            axis=2)

	nd = int(np.round(dt/delta_t, 0))
	nf = int(np.round(tfinal/delta_t, 0))
	ndac = tp.shape[0]
	# This could be a way to check that we got a list of zpk/tf/lti objects 
	# instead of a single object. Right now, we do not check. 
	#if type(S) == tuple and type(S[0]) == tuple and np.isscalar(S[0][0]):
	#	S = [[lti(S)]]
	ni = len(S) # number of inputs
	
	if ni % ndac != 0:
		raise ValueError('The number of inputs must be divisible by the number of dac timings.')
		# Original comment from the MATLAB sources:
		# This requirement comes from the complex case, where the number of inputs
		# is 2 times the number of dac timings. I think this could be tidied up.

	nis = int(ni/ndac) # Number of inputs grouped together with a common DAC timing
	                   # (2 for the complex case)
	# Notice the number of outputs is hard-coded to one in the following
	# delsigma actually allows for supplying pulse() objects that have 
	# a variable number of outputs. As we have no MIMO description available,
	# we consider only 1-output systems
	if not nosum: # Sum the responses due to each input set
		y = np.zeros((np.ceil(tfinal/float(dt)) + 1, 1, nis))
	else:
		y = np.zeros((np.ceil(tfinal/float(dt)) + 1, 1, ni))

	for i in range(ndac):
		n1 = int(np.round(tp[i, 0]/delta_t, 0))
		n2 = int(np.round(tp[i, 1]/delta_t, 0))
		z1 = (n1, y1.shape[1], nis)
		z2 = (n2, y1.shape[1], nis)
		yy = + np.concatenate((np.zeros(z1), y1[:nf-n1+1, :, i*nis:(i + 1)*nis]), axis=0) \
		     - np.concatenate((np.zeros(z2), y1[:nf-n2+1, :, i*nis:(i + 1)*nis]), axis=0)
		yy = yy[::nd, :, :]
		if not nosum: # Sum the responses due to each input set
			y = y + yy
		else:
			y[:, :, i] = yy.reshape(yy.shape[0:2])
	return y

def test_pulse():
	"""Test function for pulse()"""
	import pkg_resources
	import scipy.io
	H = ([1], [1, 2, 10])
	pp = pulse([H], tp=(0., 1.), dt=.1, tfinal=10., nosum=False)
	pp = pp.reshape((pp.shape[0],1))
	fname = pkg_resources.resource_filename(__name__, "test_data/test_pulse.mat")
	pp2 = scipy.io.loadmat(fname)['pp']
	assert np.allclose(pp, pp2, atol=1e-6, rtol=1e-4)
	# FIXME ALSO CHECK MIMO TFS

