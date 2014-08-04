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

from warnings import warn
import numpy as np
import pylab as plt

# Plot lollipops (o's and sticks)
# 
#       ^
#       | o     o     o         o o o   o o 
#       | |   o | o o | o     o | | | o | | 
#       | | o | | | | | | o o | | | | | | | 
#       +----------------------------------->

def lollipop(x, y, color=None, lw=2, ybot=0):
    """Plot lollipops (o's and sticks)
    
    **Parameters:**

    x, y : ndarrays
        The data to be plotted

    color : any matplotlib color, optional
            plotting color

    lw : float, optional
         line width value in points

    ybot : float, optional
           Dummy parameter available for compatiblity

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

    if ybot:
        warn('lollipop() got a non-zero ybot, but only ybot=0 is ' + \
             'supported. Setting ybot to 0.')
    markerline, stemlines, baseline = plt.stem(x, y, '-')
    if not color or color == 'None':
        color = stemlines[0].get_color()
    lolli_fmt = {'linewidth': lw, 'color': color}
    pop_fmt = {'mec': color, 'markerfacecolor':'None',  \
               'markersize':10, 'markeredgewidth': lw*1.1}
    plt.setp(markerline, **pop_fmt)
    plt.setp(stemlines, **lolli_fmt)
    plt.setp(baseline, 'color','k')

