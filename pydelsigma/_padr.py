# -*- coding: utf-8 -*-
# _padr.py
# This module provides the padr function.
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

"""This module provides the padr() function, which pads a matrix on the 
right.
"""

import numpy as np

def padr(x, n, val=0.):
	"""y = padb(x, n, val)
	 Pad a matrix x on the right to length n with value val
	The empty matrix is assumed to be have 1 empty row
	"""
	if len(x.shape) == 1:
		xp = x.reshape((1, x.shape[0]))
	else:
		xp = x
	y = np.concatenate(
	                   (xp, 
	                    val*np.ones((xp.shape[0], n - xp.shape[1]))
	                   ), axis=1
	                  )
	return y
	                  
def test_padr():
	"""Test function for padr()"""
	tv = np.eye(15)
	tr = padr(tv, n=25, val=2)
	res = np.concatenate((tv, 2.*np.ones((15, 10))), axis=1)
	assert np.allclose(tr, res, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_padr()
