# -*- coding: utf-8 -*-
# _partitionABCD.py
# Module providing the partitionABCD function
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

"""Module providing the partitionABCD() function
"""

def partitionABCD(ABCD, m=None):
	"""Partition ABCD into A, B, C, D for an m-input state-space system."""
	# remember the ABCD matrix is assembled like this:
	# [[A, B], 
	#  [C, D]]
	if m is None:
		n = min(ABCD.shape) - 1
		m = ABCD.shape[1] - n
	else:
		n = ABCD.shape[1] - m

	r = ABCD.shape[0] - n

	A = ABCD[:n, :n]
	B = ABCD[:n, n:n+m]
	C = ABCD[n:n+r, :n]
	D = ABCD[n:n+r, n:n+m]
	return A, B, C, D

def test_partitionABCD():
	"""Test function for partitionABCD()"""
	# FIXME: no test for m-input systems
	import numpy as np
	from scipy.signal import lti
	ob = lti((1, ), (1, 2, 10))
	ab = np.hstack((ob.A, ob.B))
	cd = np.hstack((ob.C, ob.D.reshape((1,1))))
	abcd = np.vstack((ab, cd))
	a, b, c, d = partitionABCD(abcd)
	assert np.allclose(a, ob.A, rtol=1e-5, atol=1e-8)
	assert np.allclose(b, ob.B, rtol=1e-5, atol=1e-8)
	assert np.allclose(c, ob.C, rtol=1e-5, atol=1e-8)
	assert np.allclose(d, ob.D, rtol=1e-5, atol=1e-8)
	ABCD = [[1.000000000000000, 0., 0., 0.044408783846879, -0.044408783846879],
	        [0.999036450096481, 0.997109907515262, -0.005777399147297, 0., 0.499759089304780],
	        [0.499759089304780, 0.999036450096481, 0.997109907515262,  0., -0.260002096136488],
	        [0,                 0,                 1.000000000000000,  0, -0.796730400347216]]
	ABCD = np.array(ABCD)
	at = np.array([[1., 0., 0.],
	              [0.999036450096481, 0.997109907515262, -0.005777399147297],
	              [0.499759089304780, 0.999036450096481,  0.997109907515262]
	             ])
	bt = np.array([[0.044408783846879, -0.044408783846879],
	              [0., 0.499759089304780],
                      [0., -0.260002096136488]
                     ])
	ct = np.array([0., 0, 1])
	dt = np.array([0., -0.796730400347216])
	ar, br, cr, dr = partitionABCD(ABCD)
	assert np.allclose(at, ar, rtol=1e-5, atol=1e-8)
	assert np.allclose(bt, br, rtol=1e-5, atol=1e-8)
	assert np.allclose(ct, cr, rtol=1e-5, atol=1e-8)
	assert np.allclose(dt, dr, rtol=1e-5, atol=1e-8)

