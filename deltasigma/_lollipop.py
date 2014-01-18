# -*- coding: utf-8 -*-
# _lollipop.py
# The lollipop function: plot digital samples like it was a lollipop field.
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

"""Lollipop plotting function.
"""

import numpy as np
import pylab as plt
from numpy.matlib import repmat

# Plot lollipops (o's and sticks)
# 
#	   ^
#	   | o     o     o         o o o   o o 
#	   | |   o | o o | o     o | | | o | | 
#	   | | o | | | | | | o o | | | | | | | 
#          +----------------------------------->

def lollipop(x, y, color='b', lw=2, ybot=0):
	"""Plot lollipops (o's and sticks)
	
	**Parameters:**

	x, y : ndarray
	       data to be plotted

	color : any matplotlib color
	        plotting color

	lw : float
             line width value in points

	ybot : float
               ground level for the sticks.

	**Returns:**

	None
	
	**Example:**

	.. plot::

	    import matplotlib.pyplot as plt
	    import numpy as np
	    from deltasigma import lollipop
	    t = np.arange(1, 20)*1e-3
	    f = 20.
	    a = np.sin(2*np.pi*f*t)
	    lollipop(t, a)
	    plt.gcf().set_size_inches((8, 4))
	    plt.grid(True)
	    plt.show()
	"""

	hold_status = plt.ishold()
	plt.hold(True)
	lolli_fmt = {'linewidth': lw, 'color': color}
	pop_fmt = {'markeredgecolor': color, 'markerfacecolor':'None', 'linestyle': 'None', \
	                'marker': 'o', 'markersize':10, 'markeredgewidth': lw*1.1}
	# Plot the circles
	plt.plot(x, y, **pop_fmt)

	# Make x and y row vectors, then plot as sticks
	x = x.transpose()
	y = y.transpose()
	x = np.vstack((x, x, float('NaN')*np.ones(x.shape)))
	y = np.vstack((y, repmat(ybot, 2, y.shape[0])))
	plt.plot(x, y, **lolli_fmt)

	plt.hold(hold_status)

def test_lollipop():
	"""Test function for lollipop()"""
	t = np.arange(1, 20)*1e-3
	f = 20.
	a = np.sin(2*np.pi*f*t)
	plt.figure()
	lollipop(t, a)
	plt.grid(True)
