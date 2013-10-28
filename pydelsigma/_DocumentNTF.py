# -*- coding: utf-8 -*-
# _DocumentNTF.py
# Module providing DocumentNTF
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

"""Module providing configuration switches.
"""

from __future__ import division
import numpy as np
import pylab as plt

from ._plotPZ import plotPZ
from ._ds_freq import ds_freq
from ._dbv import dbv
from ._evalTF import evalTF
from ._ds_f1f2 import ds_f1f2
from ._rmsGain import rmsGain
from ._calculateTF import calculateTF
from ._infnorm import infnorm

def DocumentNTF(arg1, osr=64, f0=0, quadrature=False):
	"""Plot the NTF's poles and zeros as well as its frequency-response

	The first argument is either the NTF or ABCD matrix. 
	If the first argument is ABCD, the stf is also plotted.
	"""
	if isinstance(arg1, np.ndarray):
		ABCD = arg1
		ntf, stf = calculateTF(ABCD)
	else:
		ntf = arg1
	
	fig = plt.figure()
	plt.subplot(211)
	plotPZ(ntf, 'b', 6)
	f = ds_freq(osr, f0, quadrature)
	z = np.exp(2j * np.pi * f)
	H = dbv(evalTF(ntf, z))
	plt.plot(f, H, 'b')
	plt.subplot(212)
	
	if 'stf' in globals() and not (0 in stf.shape):
		set_(gcf, 'name', 'NTF and STF')
		G = dbv(evalTF(stf, z))
		plt.hold(True)
		plt.plot(f, G, 'm')
		plt.hold(False)
	
	f1, f2 = ds_f1f2(osr, f0, quadrature)
	NG0 = dbv(rmsGain(ntf, f1, f2))
	plt.hold(True)
	plt.plot(np.array([f1, f2]).reshape(1, -1), NG0 * np.array([1, 1]).reshape(1, -1), 'k', linewidth=3)
	
	if f0  ==  0:
		text(0.5 / osr, NG0, sprintf('  %.0fdB', NG0), 'Vert', 'Mid', 'Hor', 'Left')
	else:
		text(f0, NG0 + 1, sprintf('%.0fdB', NG0), 'Vert', 'Bot', 'Hor', 'Cen')
	msg = ' Inf-norm of H = %.2f\n 2-norm of H = %.2f', infnorm(ntf), rmsGain(ntf, 0, 1)
	if f0 < 0.25:
		text(0.48, 0, msg, 'Vert', 'Top', 'Hor', 'Right')
	else:
		text(f_left, 0, msg, 'Vert', 'Top')
	if quadrature:
		ING0 = dbv(rmsGain(ntf, - f1, - f2))
		plt.plot(- np.array([f1, f2]).reshape(1, -1), ING0 * np.array([1, 1]).reshape(1, -1), 'k', 'Linewidth', 3)
		text(- f0, ING0 + 1, sprintf('%.0fdB', ING0), 'Vert', 'Bot', 'Hor', 'Cen')
		f_left = - 0.5
	else:
		f_left = 0
	#figureMagic(np.array([f_left, 0.5]).reshape(1, -1), 1 / 16, 2, np.array([- 80, 15]).reshape(1, -1), 10, 2)
	plt.xlabel('frequency')
	plt.title('Frequency Response')
	return fig
	
def text(*args):
	"""Dummy function for now"""
	pass
