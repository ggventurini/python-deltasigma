# -*- coding: utf-8 -*-
# _figureMagic.py
# Module providing figureMagic()
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

"""Module providing the utility function figureMagic()
"""

from __future__ import division
from warnings import warn
import numpy as np
import pylab as plt

def figureMagic(xRange=None, dx=None, xLab=None, yRange=None, dy=None, yLab=None, 
	        size=None, name=None):
	"""Utility function to quickly set plot parameters.
	Any unspecified parameter is left untouched.
	
	.. plot::

	   import numpy as np
	   import pylab as plt
	   from pydelsigma import figureMagic
	   t = np.linspace(0, 1e-3)
	   a = np.sin(2*np.pi*1e3*t + np.pi/4)
	   plt.plot(t, a)
	   figureMagic([0, 1e-3], dx=.1e-3, xLab=None, yRange=[-1.2, 1.2],
	               dy=.1, yLab=None, size=(8, 4), name="Sin wave, f = 1kHz")

	"""
	fig = plt.gcf()
	if size is not None and len(size):
		fig.set_size_inches(size)
	
	plt.grid(True)

	ax = plt.gca()
	if name is not None:
		ax.set_title(name, fontsize=14)

	if xRange is not None or yRange is not None:
		ax.set_autoscale_on(False)
		ax.set_aspect('auto', 'box')

	if xRange is not None:
		ax.set_xlim(xRange)
	if yRange is not None:
		ax.set_ylim(yRange)

	if dx is not None:
		x1, x2 = ax.get_xlim()
		ax.xaxis.set_ticks(np.arange(x1, x2, dx))
	if dy is not None:
		y1, y2 = ax.get_ylim()
		ax.yaxis.set_ticks(np.arange(y1, y2, dy))
	#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
	
	if xLab is not None:
		warn('figureMagic() got xLab=%s, but xLab is not implemented and will be ignored.' %
			 xLab)
	if yLab is not None:
		warn('figureMagic() got yLab=%s, but xLab is not implemented and will be ignored.' %
			 yLab)
	return
