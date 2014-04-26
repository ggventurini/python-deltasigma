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
    """Plot the spectrum of a band-pass modulator in dB.

    The plot is a logarithmic plot, centered in 0, corresponding to f0,
    extending to negative frequencies, with respect to the center frequencies
    and to positive frequencies.

    The plot employs a logarithmic x-axis transform far from the origin and a
    linear one close to it, allowing the x-axis to reach zero and extend to
    negative values as well.

    .. note::
        This is implemented in a slightly different way from The MATLAB Delta
        Sigma Toolbox, where all values below ``xmin`` are clipped and the scale is
        always logarithmic. It our implementation, no clippin is done and below
        ``xmin`` the data is simply plotted with a linear scale. For this reason
        slightly different plots may be generated.

    **Parameters:**

    V : 1d-ndarray or sequence
        Hann-windowed FFT

    f0 : int
        bin number of center frequency

    fbin : int
        bin number of test tone

    x : 3-elements sequence-like
        x is a sequence of three *positive* floats: ``xmin``, ``xmax_left``, ``xmax_right``.
        ``xmin`` is the minimum value of the logarithmic plot range. ``xmax_left`` is the
        length of the plotting interval on the left (negative) side, ``xmax_right`` is its
        respective on the right (positive) side.

    y : 3-elements sequence-like
        y is a sequence of three floats: ``ymin``, ``ymax``, ``dy``.
        ``ymin`` is the minimum value of the y-axis, ``ymax`` its maximum value and
        ``dy`` is the ticks spacing.

    .. note::
        The MATLAB Delta Sigma toolbox allows for a fourth option ``y_skip``, which
        is the ``incr`` value passed to MATLAB's ``axisLabels``.
        No such thing is supported here. A warning is issued if ``len(v) == 4``.

    Additional keyword parameters ``**fmt`` will be passed to matplotlib's ``semilogx()``.

    The FFT is smoothed before plotting and converted to dB. See
    :func:`logsmooth` for details regarding the algorithm used.

    **Returns:**

    *None*

    .. plot::

        from __future__ import division
        from deltasigma import synthesizeNTF, simulateDSM
        from deltasigma import calculateSNR, ds_hann, bilogplot
        import pylab as plt
        import numpy as np
        f0 = 1./8
        OSR = 64
        order = 8
        N = 8192
        H = synthesizeNTF(order, OSR, 1, 1.5, f0)
        fB = int(np.ceil(N/(2. * OSR)))
        ftest = int(np.round(f0*N + 1./3 * fB))
        u = 0.5*np.sin(2*np.pi*ftest/N*np.arange(N))
        v, xn, xmax, y = simulateDSM(u, H)
        spec = np.fft.fft(v*ds_hann(N))/(N/4)
        X = spec[:N/2 + 1]
        plt.figure()
        bilogplot(X, f0*N, ftest, (.03, .3, .3), (-140, 0, 10))

    """
    V = carray(V)
    if len(V.shape) > 1:
        if np.prod(V.shape) > max(V.shape):
            raise ValueError("The input value V should have only one" +
                             " non-unitary dimension.")
        V = V.squeeze()
    Xl = V[f0::-1]
    Xr = V[f0:]
    N = V.shape[0] - 1
    fbin = abs(fbin - f0)
    fl, pl = _logsmooth2(Xl, fbin)
    fr, pr = _logsmooth2(Xr, fbin)
    p = np.concatenate((pl[::-1], pr))
    f = np.concatenate((-fl[::-1], fr))
    plt.plot(f, p, **fmt)
    plt.xscale('symlog', linthreshx=x[0],
               subsx=np.logspace(10**int(np.ceil(np.log10(x[0]))),
                                 10**int(1+np.ceil(np.log10(max(x[2], x[1])))))
               )
    ax = plt.gca()
    ax.set_xlim([-x[1], x[2]])
    ax.set_ylim([y[0], y[1]])
    plt.grid(True)
    ytix = range(y[0], y[1] + 1, y[2])
    ax.yaxis.set_ticks(ytix)
    # we do not support axis labels
    # set_(gca,'YTickLabel', axisLabels(ytix, y[3]))
    #
    if len(y) == 4 and not y[3] is None:
        warn("Specifying y_skip is not currently supported and " +
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
    f1 = int((inBin - 1) % n) + 1
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
        p[i] = 10*np.log10(norm(X[(startbin[i] - 1):stopbin[i]])**2 /
                           (stopbin[i] - startbin[i] + 1))
    return f, p


def test_bilogplot():
    """Test function for bilogplot()"""
    from ._synthesizeNTF import synthesizeNTF
    from ._simulateDSM import simulateDSM
    from ._calculateSNR import calculateSNR
    from ._ds_hann import ds_hann
    f0 = 1./8
    OSR = 64
    order = 8
    N = 8192
    H = synthesizeNTF(order, OSR, 1, 1.5, f0)
    fB = int(np.ceil(N/(2. * OSR)))
    ftest = int(np.round(f0*N + 1./3 * fB))
    u = 0.5*np.sin(2*np.pi*ftest/N*np.arange(N))
    v, xn, xmax, y = simulateDSM(u, H)
    spec = np.fft.fft(v*ds_hann(N))/(N/4)
    X = spec[:N/2 + 1]
    plt.figure()
    bilogplot(X, f0*N, ftest, (.03, .3, .3), (-140, 0, 10))
    # graphical function: we check it doesn't fail
