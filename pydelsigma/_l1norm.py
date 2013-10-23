# -*- coding: utf-8 -*-
# _l1norm.py
# Module providing the l1norm function
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

"""Module providing the l1norm() function
"""

import numpy as np
from scipy.signal import impulse

def l1norm(H):
	"""Compute the l1-norm of a z-domain transfer function.
	"""
	_, y = impulse(H, T=np.arange(544)/5.43) # how the faq does MATLAB pick the time points?
	return np.sum(np.abs(y))

def test_l1norm():
	"""Test function for l1norm
	"""
	zeros = np.array(())
	poles = np.array((-.5,))
	k = 1.
	zpk_tuple = zeros, poles, k
	assert np.allclose(l1norm(zpk_tuple), 11.3676723351, rtol=1e-5, atol=1e-8)

