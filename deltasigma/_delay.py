# -*- coding: utf-8 -*-
# _delay.py
# This module provides the delay function.
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

"""This module provides the delay() function, used to delay a signal by a 
specified amount of samples.
"""

import numpy as np

def delay(x, n=1):
	"""Delay signal ``x`` by ``n`` samples.
	"""
	nx = np.size(x)
	if nx <= n:
		y = np.zeros(x.shape)
	else:
		if len(x.shape) == 1:
			y = np.concatenate((np.zeros((n,)), x[:nx-n]))
		elif x.shape[0] > x.shape[1]:
			y = np.concatenate((np.zeros((n, x.shape[1])), x[:nx-n, :]), axis=0)
		else:
			y = np.concatenate((np.zeros((x.shape[0], n)), x[:, :nx-n]), axis=1)
	return y

def test_delay():
	"""Test function for delay()
	"""
	v1 = np.arange(4)
	v2 = np.ones((10, 1))
	v3 = np.ones((1, 10))
	r1 = np.array((0, 0, 0, 0))
	r2 = np.zeros((10, 1))
	r2[5:, 0] = 1
	r3 = np.zeros((1, 10))
	r3[0, 5:] = 1
	assert np.allclose(delay(v1, 5), r1)
	assert np.allclose(delay(v2, 5), r2)
	assert np.allclose(delay(v3, 5), r3)

