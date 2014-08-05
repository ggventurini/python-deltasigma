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

"""Module providing the DocumentNTF function.
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
from ._figureMagic import figureMagic

def DocumentNTF(arg1, osr=64, f0=0, quadrature=False):
    """Plot the NTF's poles and zeros as well as its frequency-response

    The first argument is either the NTF or ABCD matrix.
    If the first argument is ABCD, the STF is also plotted.

    .. plot::

       from deltasigma import *
       import numpy as np
       import pylab as plt
       order = 4
       osr = 64
       nlev = 2
       f0 = 0.
       Hinf = 1.5
       form = 'CRFB'
       ntf = synthesizeNTF(order, osr, 2, Hinf, f0)
       a, g, b, c = realizeNTF(ntf, form)
       b = np.concatenate(( # Use a single feed-in for the input
                           np.atleast_1d(b[0]),
                           np.zeros((b.shape[0] - 1,))
                         ))
       ABCD = stuffABCD(a, g, b, c, form)
       DocumentNTF(ABCD, osr, f0)

    """
    if isinstance(arg1, np.ndarray):
        ABCD = arg1
        ntf, stf = calculateTF(ABCD)
    else:
        ntf = arg1
        stf = None

    fig = plt.figure(figsize=(12, 5))
    plt.subplot(121)
    plotPZ(ntf, 'b', 6, showlist=False)
    plt.title('Poles and Zeros')
    plt.subplot(122)
    f = ds_freq(osr, f0, quadrature)
    z = np.exp(2j * np.pi * f)
    H = dbv(evalTF(ntf, z))
    plt.plot(f, H, 'b')

    if stf is not None:
        fig.suptitle('NTF and STF', fontsize=14)
        G = dbv(evalTF(stf, z))
        plt.hold(True)
        plt.plot(f, G, 'm')
        plt.hold(False)
    else:
        fig.suptitle('NTF', fontsize=14)

    f1, f2 = ds_f1f2(osr, f0, quadrature)
    NG0 = dbv(rmsGain(ntf, f1, f2))
    plt.hold(True)
    plt.plot(np.array([f1, f2]), NG0*np.array([1, 1]), 'k', linewidth=3)

    if f0  ==  0:
        plt.text(0.5/osr, NG0, '  %.0fdB' % NG0, horizontalalignment = 'left',
                 verticalalignment = 'center')
    else:
        plt.text(f0, NG0 + 1, '%.0fdB' % NG0, horizontalalignment = 'center',
                 verticalalignment = 'bottom')
    msg = ' Inf-norm of H = %.2f\n 2-norm of H = %.2f' % (infnorm(ntf)[0], rmsGain(ntf, 0, 1))
    if quadrature:
        ING0 = dbv(rmsGain(ntf, - f1, - f2))
        plt.plot(-np.array([f1, f2]), ING0*np.array([1, 1]), 'k', linewidth=3)
        plt.text(-f0, ING0 + 1, '%.0fdB' % ING0, horizontalalignment = 'center',
                 verticalalignment = 'bottom')
        f_left = - 0.5
    else:
        f_left = 0
    # variable 'f_left' used before assignment in DocumentNTF.m #REP
    if f0 < 0.25:
        plt.text(0.48, 0, msg, horizontalalignment = 'right', verticalalignment = 'top')
    else:
        plt.text(f_left, 0, msg, horizontalalignment = 'left', verticalalignment = 'top')
    plt.grid(True)
    y_bot = min(-80, np.round(NG0*1.1, -1))
    figureMagic(xRange=(f_left, 0.5), dx=1/20., yRange=(y_bot, 15), dy=10)
    plt.ylabel('|H(f)| dB')
    plt.xlabel('Frequency ($1 \\rightarrow f_{s}$)')
    plt.title('Frequency Response')
    return

