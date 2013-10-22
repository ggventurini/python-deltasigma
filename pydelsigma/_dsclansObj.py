# -*- coding: utf-8 -*-
# _dsclansObj.py
# Module providing the dsclansObj function
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

"""Module providing the dsclansObj() function
"""

from warnings import warn
import numpy as np
from scipy.signal import impulse
from ._evalTF import evalTF
from ._dsclansNTF import dsclansNTF

def dsclansObj(x, order, OSR, Q, rmax, Hz):
	"""Objective function for clans.m 
	f is the magnitude of H at the band-edge
	g = ||h||_1 - Q
	"""

	# Translate x into H.
	H = dsclansNTF(x, order, rmax, Hz)
	warn("Warning untested function.")
	# Compute f and g
	f = np.abs(evalTF(H, np.exp(np.pi*1.j/OSR)))
	T = np.arange(101)
	y, _ = impulse((H.zeros, H.poles, H.k), T=T)
	g = np.sum(np.abs(y)) - 1. - Q

	return f, g

def test_dsclansObj():
	"""Test function for dsclansNTF()
	"""
	OSR = 64.
	Q = 10.
	x = np.ones((10,))
	order = 3
	rmax = 5.
	Hz = 100.
	dsclansObj(x, order, OSR, Q, rmax, Hz)
	assert False

	
	
