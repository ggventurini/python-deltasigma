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
import pylab as plt

from ._utils import _get_zpk

def plotPZ(H, color='b', markersize=5, showlist=False):
    """Plot the poles and zeros of a transfer function.

    **Parameters:**

    H : transfer function
        Any supported transfer function representation, 
        eg num/den, zpk, lti...

    color : Any matplotlib-compatible color descr, optional
        For example, 'r' for 'red' or '#000080' for 'navy'.
        You can also specify separately poles and zeros, in a tuple.

    markersize : scalar, optional
        The markers size in points.

    showlist : boolean, optional
        Superimpose a list of the poles and zeros on the plot.

    .. plot::

       import pylab as plt
       from deltasigma import synthesizeNTF, plotPZ
       order = 5
       osr = 32
       f0 = 0.
       Hinf = 1.5
       ntf = synthesizeNTF(order, osr, 2, Hinf, f0)
       plt.figure(figsize=(8, 6))
       plotPZ(ntf, color=('r', 'b'), showlist=True)
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
    
    z, p, _ = _get_zpk(H)

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
    
    # plt.axes().set_aspect('equal', 'datalim')
    plt.ylabel('Imag')
    plt.xlabel('Real')

    if not hold_status:
        plt.hold(False)

