# -*- coding: utf-8 -*-
# _plotPZ.py
# Module providing the plotPZ function
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

"""Module providing plotPZ(), useful to plot the poles and zeros of a transfer function.
"""

import numpy as np
from scipy.signal import tf2zpk
import pylab as plt

def plotPZ(H, color='b', markersize=5, showlist=False):
	"""function plotPZ(H,color='b',markersize=5,list=0)
	Plot the poles and zeros of a transfer function.
	If list is non-zero, a list of the poles and zeros is superimposed on the plot.
	"""

	# Parts of the code come from 'pydsm'
	#
	# Original copyright notices:
	#
	# For pydsm
	# Copyright (c) 2012, Sergio Callegari
	# All rights reserved.
	#
	# For Richard Schreier's Delta Sigma toolbox
	# Copyright (c) 2009, Richard Schreier
	
	if (hasattr(H, 'inputs') and not H.inputs == 1) or \
	   (hasattr(H, 'outputs') and not H.outputs == 1):
			raise TypeError, "Only SISO transfer functions can be evaluated."
	if hasattr(H, 'num') and hasattr(H, 'den'):
		filt = hasattr(H, 'outputs')
		num = H.num[0][0] if filt else H.num
		den = H.den[0][0] if filt else H.den
		z, p, k = tf2zpk(num, den)
	elif (hasattr(H, 'zeros') and hasattr(H, 'poles')) or \
	   (hasattr(H, 'zero') and hasattr(H, 'pole')):
		# LTI objects have poles and zeros, 
		# TransferFunction-s have pole() and zero()
	   	z = H.zeros if hasattr(H, 'zeros') else H.zero()
	   	p = H.poles if hasattr(H, 'poles') else H.pole()
		if hasattr(H, 'k'):
			k = H.k
		elif hasattr(H, 'gain'):
			k = H.gain  
		elif hasattr(H, 'returnScipySignalLti'): 
			k = np.array(H.returnScipySignalLti()[0][0].gain)
	elif hasattr(H, 'form') and H.form == 'zp':
		z, p, k = H.k, H.zeros, H.poles
	elif hasattr(H, 'form') and H.form == 'coeff':
		z, p, k = tf2zpk(H.num, H.den)
	elif hasattr(H, 'form'):
		raise ValueError, '%s: Unknown form: %s' % (__name__, H.form)
	elif hasattr(H, '__len__'):
		if len(H) == 2:
			z, p, k = tf2zpk(H[0], H[1])
		elif len(H) == 3:
			z, p, k = H
	else:
		raise TypeError, '%s: Unknown transfer function %s' % (__name__, str(H))

	pole_fmt = {'marker': 'x', 'markersize': markersize}
	zero_fmt = {'marker': 'o', 'markersize': markersize}

	if isinstance(color, list) or isinstance(color, tuple):
		pole_fmt['color'] = color[0]
		zero_fmt['color'] = color[1]
	else:
		pole_fmt['color'] = color
		zero_fmt['color'] = color

	hold_status = plt.ishold()
	plt.grid(True)

	# Plot x and o for poles and zeros, respectively
	plt.plot(p.real, p.imag, linestyle='None', **pole_fmt)
	plt.hold(True)
	if len(z) > 0:
		plt.plot(z.real, z.imag, linestyle='None', **zero_fmt)

	# Draw unit circle, real axis and imag axis
	circle = np.exp(2j*np.pi*np.linspace(0, 1, 100))
	plt.plot(circle.real, circle.imag)
	plt.axis('equal')

	if showlist:
		# List the poles and zeros
		pp = p[p.imag >= 0]
		y = 0.05*(len(pp)+1)
		str_p = 'Poles:'
		plt.text(-0.9, y, str_p,
			 horizontalalignment = 'left',
			 verticalalignment = 'center')
	y = y - 0.1
	for i in xrange(len(pp)):
	    if pp[i].imag == 0:
		str_p = '%+.4f' % pp[i].real
	    else:
		str_p = '%+.4f+/-j%.4f' %  (pp[i].real, pp[i].imag)
	    plt.text(-0.9, y, str_p,
		     horizontalalignment = 'left',
		     verticalalignment = 'center')
	    y = y - 0.1
	if len(z) > 0:
	    zz = z[z.imag >= 0]
	    y = 0.05*(len(zz)+1)
	    str_z = 'Zeros:'
	    plt.text(0, y, str_z,
		     horizontalalignment = 'left',
		     verticalalignment = 'center')
	    y = y - 0.1
	    for i in xrange(len(zz)):
		if zz[i].imag == 0:
		    str_z = '%+.4f' % zz[i].real
		else:
		    str_z = '%+.4f+/-j%.4f' % (zz[i].real, zz[i].imag)
		plt.text(0, y, str_z,
		         horizontalalignment = 'left',
		         verticalalignment = 'center')
		y = y - 0.1

	plt.ylabel('Imag')
	plt.xlabel('Real')

	if not hold_status:
		plt.hold(False)

if __name__ == '__main__':
	plt.figure()
	plotPZ(((1, .2), (1, 2, .10)), showlist=True)
	plt.show()
