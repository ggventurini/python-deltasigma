# -*- coding: utf-8 -*-
# dbm.py
# This module provides the dbm function.
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

from __future__ import division
import numpy as np

def dbm(v, R=50):
	""" dbm(v, R=50) = 10*log10(v^2/50*1000)  
	The equivalent in dBm of an rms voltage v
	"""
	if not hasattr(v, 'shape'):
		if not hasattr(v, '__len__'):
			v = np.array((v,))
		else:
			v = np.array(v)
	elif v.shape == ():
		v = np.array((v,))
	y = -np.Inf*np.ones(np.size(v))
	nonzero = (v != 0)
	y[nonzero] = 10.*np.log10(np.abs(v[nonzero]**2.)/R) + 30
	return y

def test_dbm():
	v = np.arange(10)*1e-3
	r = [-np.inf, -46.98970004, -40.96910013, -37.44727495, -34.94850022,
	     -33.01029996, -31.42667504, -30.08773924, -28.9279003, -27.90484985]
	assert np.allclose(dbm(v), r, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_dbm()
