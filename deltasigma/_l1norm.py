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

from __future__ import division
import numpy as np
from scipy.signal import dimpulse

def l1norm(H):
	"""Compute the l1-norm of a z-domain transfer function.
	"""
	_, y = dimpulse(H, t=np.arange(100))
	return np.sum(np.abs(y[0]))

def test_l1norm():
	"""Test function for l1norm()
	"""
	zeros = np.array(())
	poles = np.array((.5,))
	k = 1.
	zpkt_tuple = zeros, poles, k, 1
	assert np.allclose(l1norm(zpkt_tuple), 2., rtol=1e-5, atol=1e-8)

