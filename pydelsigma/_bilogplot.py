# -*- coding: utf-8 -*-
# _bilogplot.py
# Module providing the bilogplot function
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

"""Module providing the bilogplot() function
"""

from __future__ import division

import numpy as np
import pylab as plt

from warnings import warn
from numpy.linalg import norm

from ._dbp import dbp
from ._utils import carray

def bilogplot(V, f0, fbin, x, y, **fmt):
    """Plot two side-by-side spectra

    **Parameters:**

    V : 1d-ndarray
        Hann-windowed FFT

    f0 : int
        bin number of center frequency

    fbin : int
        bin number of test tone

    x : 3-elements sequence-like
        x is a sequence of `xmin, xmax_left, xmax_right`.

    y : 3-elements sequence-like
        y is a sequence of `ymin, ymax, dy, y_skip`.

    Additional keyword parameters `**fmt` will be passed to matplotlib's semilogx()
    
    The FFT is smoothed before plotting and converted to dB. See
    :func:`logsmooth` for details regarding the algorithm used. 

    **Returns:**

    *None*

    """
    V = carray(V)
    if len(V.shape) > 1:
        if np.prod(V.shape) > max(V.shape):
            raise ValueError("The input value V should have only one" + \
                            " non-unitary dimension.") 
        V = V.squeeze()
    Xl = V[f0:-1]
    Xr = V[f0:]
    N = V.shape[0] - 1
    fbin = abs(fbin - f0)
    f, p = _logsmooth2(Xl, fbin)
    f = f0/N*f
    # SUBPLOT #1
    plt.subplot(211)
    plt.semilogx(f, p, **fmt)
    plt.hold(True)
    ax = plt.gca()
    # reverse the x-axis direction 
    ax.xlim([x[0], x[1]])
    ax.invert_xaxis()
    ax.ylim([y[0], y[1]])
    plt.grid(True)
    ytix = range(y[0], y[1] + 1, y[2])
    ax.yaxis.set_ticks(ytix)
    # we do not support axis labels
    #set_(gca,'YTickLabel', axisLabels(ytix, y[3]))
    # SUBPLOT #2
    f, p = _logsmooth2(Xr, fbin)
    f = (N - f0)/N*f
    plt.subplot(212)
    plt.semilogx(f, p, **fmt)
    plt.hold(True)
    ax = plt.gca()
    ax.xlim([x[0], x[2]])
    ax.ylim([y[0], y[1]])
    plt.grid(True)
    ytix = range(y[0], y[1] + 1, y[2])
    ax.yaxis.set_ticks(ytix)
    ax.yaxis.set_label_position("right")
    if len(y) == 4 and not y[3] is None:
        warn("Specifying ytick labels is not currently supported and " + \
             "it will be ignored. Sorry!")
    return

def _logsmooth2(X, inBin, nbin=8):
    """Smooth the fft, X, and convert it to dB.
    Use ``nbin`` bins from 0 to 3*inBin, 
    thereafter increase bin sizes by a factor of 1.1, staying less than 2^10.
    For the 3 sets of bins inBin+[0:2], 2*inBin+[0:2], and 
    3*inBin+[0:2], don't do averaging. This way, the noise BW
    and the scaling of the tone and its harmonics are unchanged.
    Unfortunately, harmonics above the third appear smaller than they 
    really are because their energy is averaged over several bins.
    """    
    N = max(X.shape)
    n = nbin
    f1 = int((inBin - 1) %  n) + 1
    startbin = np.concatenate((np.arange(f1, inBin, n),
                               np.arange(inBin, inBin + 3),
                               np.arange(inBin + 3, 2*inBin, n),
                               2*inBin + np.arange(0, 3),
                               np.arange(2*inBin + 3, 3*inBin, n),
                               3*inBin + np.arange(0, 3)))
    m = startbin[-1] + n
    while m < N:
        startbin = np.concatenate((startbin, np.array((m,))))
        n = min(n*1.1, 2**10)
        m = np.round(m + n)
    stopbin = np.concatenate((startbin[1:] - 1, np.array((N,))))
    f = ((startbin + stopbin)/2. - 1)/N
    p = np.zeros(f.shape)
    for i in range(max(f.shape)):
        p[i] = dbp(norm(X[(startbin[i] - 1):stopbin[i]])**2 /
                   (stopbin[i] - startbin[i] + 1))
    return f, p
