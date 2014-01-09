# -*- coding: utf-8 -*-
# _nabsH.py
# This module provides the nabsH function.
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

"""This module provides the nabsH() function, which computes the negative of 
the absolute value of H(z).
"""

import numpy as np

from ._evalTF import evalTF

def nabsH(w, H):
	"""Computes the negative of the absolute value of H 
	at the specified angular frequency w on the unit circle.

	This function is used by :func:`infnorm`.
	"""
	z = np.exp(1j*w)
	return -np.abs(evalTF(H, z))
	
def test_nabsH():
	"""Test function for nabsH()"""
	H = ([1, 2], [2, 0, .25], 1)
	N = 129
	w = np.linspace(0, 2*np.pi, num=N, endpoint=True)
	z = np.exp(1j*w)
	r1 = -np.abs(evalTF(H, z))
	r2 = nabsH(w, H)
	assert np.allclose(r1, r2, atol=1e-8, rtol=1e-5)

