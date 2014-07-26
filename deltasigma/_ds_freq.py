# -*- coding: utf-8 -*-
# _ds_freq.py
# This module provides the ds_freq function.
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

"""This module provides the ds_freq() function, used to generate a frequency 
vector suitable for plotting the frequency response of an NTF.
"""

from __future__ import division
import numpy as np

def ds_freq(osr=64., f0=0., quadrature=False):
	"""Frequency vector suitable for plotting the frequency response of an NTF
	"""
	if quadrature:
		f_left = -0.5
		f_special = (f0, -f0)
	else:
		f_left = 0.
		f_special = (f0, )
	f = np.linspace(f_left, 0.5, num=100)
	# Use finer spacing in the vicinity of the passband
	for fx in f_special:
		f1 = max(f_left, fx - 1./osr)
		f2 = min(0.5, fx + 2./osr)
		dels = np.where(np.logical_and(f <= f2, f >= f1))
		f = np.delete(f, dels)
		f = np.sort(np.concatenate((f, np.linspace(f1, f2, num=100))))
	return f

