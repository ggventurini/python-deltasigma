# -*- coding: utf-8 -*-
# _sinc_decimate.py
# This module provides the sinc_decimate function.
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

"""This module provides the sinc_decimate() function, which decimates a vector 
by a sinc filter of specified order and length.
"""

from __future__ import division
import numpy as np

def sinc_decimate(x, m, r):
	"""Decimate ``x`` by an ``m``-th order sinc filter of length ``r``.
	"""
	x = x[:]
	for i in range(m):
		x = np.cumsum(x)
		x = np.concatenate((x[:r], x[r:] - x[:-r]), axis=0)/r
	return x[r::r]

def test_sinc_decimate():
	"""Test function for sinc_decimate()"""
	tv = np.sin(2*np.pi*.1*np.arange(0, 100))
	r = [-1.33907057e-17, -3.55951662e-17, -2.28847549e-18, 
	     -3.55951662e-17, -1.02208548e-16, -4.35275455e-16, 
	      4.97311886e-16, -4.57479916e-16, 4.30698504e-16]
	assert np.allclose(sinc_decimate(tv, 1, 10), r, rtol=1e-5, atol=1e-8)

