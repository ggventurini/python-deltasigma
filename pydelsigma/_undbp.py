# -*- coding: utf-8 -*-
# undbp.py
# This module provides the undbp function.
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

"""This module provides the undbp() function.
"""

import numpy as np

def undbp(x):
	""" y = undbp(x)	
	Convert x from dB to a power"""
	return 10.**(x/10.)
	
def test():
	assert np.allclose([undbp(53.05)], [201836.636368], rtol=1e-05, atol=1e-08)
	assert np.allclose([undbp(3)], [1.99526231497], rtol=1e-05, atol=1e-08)

if __name__ == '__main__':
	test()
