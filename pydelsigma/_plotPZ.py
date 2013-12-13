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
	"""function plotPZ(H, color='b',markersize=5, list=0)
	Plot the poles and zeros of a transfer function.
	If list is non-zero, a list of the poles and zeros is superimposed on the plot.

	.. plot::

	   import pylab as plt
	   from pydelsigma import synthesizeNTF, plotPZ
	   order = 5
	   osr = 32
	   f0 = 0.
	   Hinf = 1.5
	   ntf = synthesizeNTF(order, osr, 2, Hinf, f0)
	   plt.figure(figsize=(8, 6))
	   plotPZ(ntf, showlist=True)
	   plt.title("NTF singularities")
	   plt.show()

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
		# for now we support both TransferFunction objects (python-control)
		# and lti objects (scipy).
		filt = hasattr(H, '__class__') and H.__class__.__name__ == 'TransferFunction'
		num = H.num[0][0] if filt else H.num
		den = H.den[0][0] if filt else H.den
		z, p, _ = tf2zpk(num, den)
	elif (hasattr(H, 'zeros') and hasattr(H, 'poles')) or \
	   (hasattr(H, 'zero') and hasattr(H, 'pole')):
		# LTI objects have poles and zeros, 
		# TransferFunction-s have pole() and zero()
		z = H.zeros if hasattr(H, 'zeros') else H.zero()
		p = H.poles if hasattr(H, 'poles') else H.pole()
	elif hasattr(H, 'form') and H.form == 'zp':
		z, p, _ = H.k, H.zeros, H.poles
	elif hasattr(H, 'form') and H.form == 'coeff':
		z, p, _ = tf2zpk(H.num, H.den)
	elif hasattr(H, 'form'):
		raise ValueError, '%s: Unknown form: %s' % (__name__, H.form)
	elif hasattr(H, '__len__'):
		if len(H) == 2:
			z, p, _ = tf2zpk(H[0], H[1])
		elif len(H) == 3:
			z, p = H[0:2]
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
	
	ax = plt.gca()
	ax.set_autoscale_on(False)
	if showlist:
		ax = plt.gca()
		x1, x2, y1, y2 = ax.axis()
		x2 = np.round((x2 - x1)*1.48 + x1, 1)
		ax.axis((x1, x2, y1, y2))
		markers = [] 
		descr = []
		ps = p[p.imag >= 0]
		for pi in ps:
			markers += [plt.Line2D((), (), linestyle='None', **pole_fmt)]
			if pi.imag == 0:
				descr += ['%+.4f' % pi.real]
			else:
				descr += ['%+.4f+/-j%.4f' %  (pi.real, pi.imag)]
		if len(z) > 0:
			for zi in z[z.imag >= 0]:
				markers += [plt.Line2D((), (), linestyle='None', **zero_fmt)]
				if zi.imag == 0:
					descr += ['%+.4f' % zi.real]
				else:
					descr += ['%+.4f +/-j%.4f' % (zi.real, zi.imag)]
		plt.legend(markers, descr, title="Poles (x) and zeros (o)", ncol=1, loc='best', 
		           handlelength=.55, prop={'size':10})
	else:
		plt.xlim((-1.1, 1.1))
		plt.ylim((-1.1, 1.1))
	plt.gca().set_aspect('equal')
	
	#plt.axes().set_aspect('equal', 'datalim')
	plt.ylabel('Imag')
	plt.xlabel('Real')

	if not hold_status:
		plt.hold(False)

if __name__ == '__main__':
	plt.figure()
	plotPZ(((1, .2), (1, 0, .10)), showlist=True)
	plt.figure()
	plotPZ(((1, .2), (1, 0, .10)), showlist=False)
	plt.show()
