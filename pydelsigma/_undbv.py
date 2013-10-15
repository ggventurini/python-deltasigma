# -*- coding: utf-8 -*-
# undbv.py
# This module provides the undbv function.
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

"""This module provides the undbv() function.
"""

import numpy as np

def undbv(x):
	""" y = undbv(x)	
	Convert x from dB to a voltage
	"""
	return 10.**(x/20.)

	
def test():
	assert np.allclose([undbv(53.05)], [449.26232467], rtol=1e-05, atol=1e-08)
	assert np.allclose([undbv(3)], [1.41253754462], rtol=1e-05, atol=1e-08)

if __name__ == '__main__':
	test()
