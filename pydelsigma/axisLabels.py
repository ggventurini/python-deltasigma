# -*- coding: utf-8 -*-
# axisLabels.py
# Module providing the axisLabel function
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

import numpy as np

def axisLabels(ran, incr):
	"""function s = axisLabels(ran, incr)
	range is an array containing the axis points (floats)
	incr may be:
	* an int, the function returns an array of strings corresponding to:
	each element of range[0]:range[-1]:incr formatted as '%g'
	* a list, the function returns an array of strings corresponding to:
	each element of incr[1]:range[-1]:incr[0] formatted as '%g'
	
	Note: all elements in ran less than 1e-6 are rounded down to 0.
	"""
	ran[np.abs(ran)<1e-6] = 0
	s = []
	if not hasattr(incr, '__len__'):
		incr = int(incr)
		first = 0
	elif len(incr) == 2:
		first = incr[1]
		incr = incr[0]
	else:
		raise ValueError, "Unrecognised incr: "+str(incr)
	for i in range(first, len(ran), incr):
		s += ['%g' % ran[i]]
	return s
	
def test_axisLabels():
	ran = np.arange(100)
	ss = axisLabels(ran, incr=10)
	r = ['0', '10', '20', '30', '40', '50', '60', '70', '80', '90']
	assert r == ss
	
if __name__ == '__main__':
	test_axisLabels()
