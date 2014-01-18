# -*- coding: utf-8 -*-
# _rmsGain.py
# Module providing the rmsGain function
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

"""Module providing the rmsGain() function
"""

import numpy as np
from scipy.linalg import norm
from ._evalTF import evalTF

def rmsGain(H, f1, f2, N=100):
	"""Compute the root mean-square gain of the discrete-time
	TF ``H`` in the frequency band ``(f1, f2)``.
	"""

	w = np.linspace(2*np.pi*f1, 2*np.pi*f2, N)
	g = norm(evalTF(H, np.exp(1j*w))) / np.sqrt(N)

	return g
	
def test_rmsGain():
	"""Test function for rmsGain()
	"""
	from ._utils import empty
	H = empty()
	H.num = (1,)
	H.den = (1, 2, 10)
	f1 = 0.001
	f2 = 0.5
	res = rmsGain(H, f1, f2, N=1000)
	res1 = 0.102245275091
	assert np.allclose((res,), (res1,), rtol=1e-05, atol=1e-08)

