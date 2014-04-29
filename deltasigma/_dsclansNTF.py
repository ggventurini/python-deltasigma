# -*- coding: utf-8 -*-
# _dsclansNTF.py
# Module providing the dsclansNTF function
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

"""Module providing the dsclansNTF() function
"""

import numpy as np
from ._utils import carray

def dsclansNTF(x, order, rmax, Hz):
	""" Conversion of clans parameters into a NTF.

	Translate x into H.
	I've changed the relationships between (zeta, wn) and x
	in order to guarantee LHP roots of the s-polynomial.
	
	Returns the NTF, a zpk tuple.
	"""
	x = x.squeeze()
	Hz = carray(Hz)
	Hz = Hz.reshape((-1,))
	Hp = np.zeros((1,), dtype=np.complex128)
	odd = (order % 2 == 1)
	if odd:
		s = -x[0]**2.
		Hp[0] = rmax*(1. + s)/(1. - s)

	for i in range(0+1*odd, order, 2):
		Hp = np.hstack((Hp, np.zeros((2,))))
		zeta = x[i]**2
		wn = x[i + 1]**2
		s = np.roots(np.array((1, 2*zeta*wn, wn**2)))
		Hp[i:i+2] = rmax*(1. + s)/(1. - s)

	H = (Hz, Hp, 1.)
	return H

def test_dsclansNTF():
	"""Test function for dsclansNTF()
	"""
	x = dsclansNTF(np.arange(1, 100.001, .001), 3, .5, 100)
	rt = np.array(x[1][:])
	rt.sort()
	rp = np.array([0, -0.016805373426715, 0.014809370415763])
	rp.sort()
	rz = (100., )
	assert np.allclose(rp, rt, rtol=1e-5, atol=1e-6)
	assert np.allclose(rz, x[0], rtol=1e-5, atol=1e-8)
	rp = np.array([0.35378443+0.j,  0.34187718+0.22831019j,
                       0.34187718-0.22831019j,  0.32978826+0.59355161j,
                       0.32978826-0.59355161j])
	x = np.array([0.67623674,  0.9277613 ,  0.70365961,  0.60374009,  0.78008118])
	Hz = np.array([0.99604531 - 0.08884669j,  0.99604531 +0.08884669j,
                       0.99860302 - 0.05283948j,  0.99860302 +0.05283948j,  1.00000000 +0.j])
	tp =  dsclansNTF(x, order=5, rmax=.95, Hz=Hz)[1]
	assert np.allclose(rp, tp, rtol=1e-5, atol=1e-6)
	
