# -*- coding: utf-8 -*-
# _lollipop.py
# The lollipop function: plot digital samples like it was a lollipop field.
# that do not find a direct replacement in numpy/scipy.
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

"""Lollipop plotting functions.
"""

import numpy as np
import pylab as plt
from numpy.matlib import repmat

def lollipop(x, y, color='b', lw=2, ybot=0):
	"""Plot lollipops (o's and sticks)
	
	o     o     o         o o o   o o
	|   o | o o | o     o | | | o | |
	| o | | | | | | o o | | | | | | |
	
	Returns:
	========
	None
	
	"""

	hold_status = plt.ishold()
	plt.hold(True)
	lolli_fmt = {'linewidth': lw, 'color': color}
	pop_fmt = {'markeredgecolor': color, 'markerfacecolor':'None', 'linestyle': 'None', \
	                'marker': 'o', 'markersize':10, 'markeredgewidth': lw*1.1}
	# Plot circles
	plt.plot(x, y, **pop_fmt)

	# Make x and y row vectors, then plot as sticks
	x = x.transpose()
	y = y.transpose()
	x = np.vstack((x, x, float('NaN')*np.ones(x.shape)))
	y = np.vstack((y, repmat(ybot, 2, y.shape[0])))
	plt.plot(x, y, **lolli_fmt)

	plt.hold(hold_status)

if __name__ == '__main__':
	t = np.arange(10)
	lollipop(t, t)
	plt.show()
