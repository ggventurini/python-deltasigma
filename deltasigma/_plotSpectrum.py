# -*- coding: utf-8 -*-
# _plotSpectrum.py
# Module providing the plotSpectrum function
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

"""Module providing the plotSpectrum() function
"""
from __future__ import division
import pylab as plt

from ._logsmooth import logsmooth

def plotSpectrum(X, fin, fmt='-', **xargs):
    """Plot a smoothed spectrum on a LOG x-axis.

    **Parameters:**

    X : 1D ndarray
        The FFT to be smoothed and plotted, *dual sided*.

    fin : int
        The bin corresponding to the input sine wave.

    fmt : string, optional
        Formatting to be passed to matplotlib's ``semilogx()``.

    **xargs : dict, optional
        Extra arguments to be passed to matplotlib's ``semilogx()``.

    Plotting is performed on the current figure.

    .. seealso::

         * :func:`logsmooth` for more information on the smoothing algorithm

         * :func:`circ_smooth` for a smoothing algorithm suitable for linear x-axes.

    .. plot::

        import numpy as np
        from deltasigma import *
        from numpy.fft import fft
        f0 = 0
        osr = 32
        quadrature = False
        Hinf = 1.5
        order = 3
        ntf = synthesizeNTF(order, osr, 0, Hinf, f0)
        f1, f2 = ds_f1f2(osr, f0, quadrature)
        delta = 2
        Amp = undbv(-3)
        f = 0.3
        N = 2**12
        f1_bin = np.round(f1*N)
        f2_bin = np.round(f2*N)
        fin = np.round(((1 - f)/2*f1 + (f + 1)/2*f2) * N)
        t = np.arange(0, N)
        u = Amp*np.cos((2*np.pi/N)*fin*t)
        v, xn, xmax, y = simulateDSM(u, ntf, 2)
        window = ds_hann(N)
        NBW = 1.5/N
        spec0 = fft(v * window)/(N/4)
        freq = np.linspace(0, 0.5, N/2 + 1)
        plt.plot(freq, dbv(spec0[:N/2 + 1]), 'c', linewidth=1, label='$S$')
        plt.hold(True)
        plotSpectrum(spec0, fin, 'b', linewidth=2.5, label='$\\\\mathrm{plotSpectrum}(S)$')
        Snn = np.abs(evalTF(ntf, np.exp(2j*np.pi*freq)))**2 * 2/12*(delta)**2
        plt.plot(freq, dbp(Snn*NBW), 'm', linewidth=1.5, label='$\\\\mathrm{from\\\\ NTF}$')
        snr = calculateSNR(spec0[f1_bin:f2_bin + 1], fin - f1_bin)
        msg = 'SQNR  =  %.1fdB\\n @ A = %.1fdBFS & osr = %.0f\\n' % \\
              (snr, dbv(spec0[fin]), osr)
        plt.text(f0 + .45/osr, -20, msg, horizontalalignment='left',
                 verticalalignment='center')
        plt.text(0.5, -3, 'NBW = %.1e ' % NBW, horizontalalignment='right',
                 verticalalignment='top')
        figureMagic((1e-3, 0.5), None, None, (-140, 0), 20, None)
        plt.ylabel('Spectrum [dB]')
        ax = plt.gca()
        ax.set_title('Third order modulator')
        plt.xlabel('Normalized frequency ($f_s \\\\rightarrow 1$)')
        plt.legend(loc=4)

    """
    f, p = logsmooth(X, fin)
    plt.semilogx(f, p, fmt, **xargs)

def test_plotSpectrum():
    """Test function for plotSpectrum()"""
    import numpy as np
    from ._ds_f1f2 import ds_f1f2
    from ._undbv import undbv
    from ._dbv import dbv
    from ._dbp import dbp
    from ._synthesizeNTF import synthesizeNTF
    from ._simulateDSM import simulateDSM
    from ._figureMagic import figureMagic
    from ._ds_hann import ds_hann
    from ._circ_smooth import circ_smooth
    from ._evalTF import evalTF
    from numpy.fft import fft
    f0 = 0
    osr = 32
    quadrature = False
    Hinf = 1.5
    order = 3
    ntf = synthesizeNTF(order, osr, 0, Hinf, f0)
    f1, f2 = ds_f1f2(osr, f0, quadrature)
    delta = 2
    Amp = undbv(-3)
    f = 0.3
    N = 2**12
    f1_bin = np.round(f1*N)
    f2_bin = np.round(f2*N)
    fin = np.round(((1 - f)/2*f1 + (f + 1)/2*f2) * N)
    t = np.arange(0, N)
    u = Amp*np.cos((2*np.pi/N)*fin*t)
    v, xn, xmax, y = simulateDSM(u, ntf, 2)
    window = ds_hann(N)
    NBW = 1.5/N
    spec0 = fft(v * window)/(N/4)
    freq = np.linspace(0, 0.5, N/2 + 1)
    # plotting
    plt.subplot(211)
    plt.plot(freq, dbv(spec0[:N/2 + 1]), 'c', linewidth=1, label='$S$')
    plt.hold(True)
    spec_smoothed = circ_smooth(np.abs(spec0)**2., 16)
    plt.plot(freq, dbp(spec_smoothed[:N/2 + 1]), 'b--', linewidth=2, label='$\\mathrm{circ\\_smooth}(S)$')
    plotSpectrum(spec0, fin, 'r', linewidth=2, label='$\\mathrm{plotSpectrum}(S)$')
    Snn = np.abs(evalTF(ntf, np.exp(2j*np.pi*freq)))**2 * 2/12*(delta)**2
    plt.plot(freq, dbp(Snn*NBW), 'm', linewidth=1.5, label='$\mathrm{from\\ NTF}$')
    plt.text(0.5, -3, 'NBW = %.1e ' % NBW, horizontalalignment='right',
             verticalalignment='top')
    figureMagic((0, 0.5), None, None, (-140, 0), 20, None)
    plt.ylabel('Spectrum [dB]')
    ax = plt.gca()
    ax.set_title('Smoothing and plotting for LOG and LIN axes')
    plt.legend(loc=4)
    plt.subplot(212)
    plt.plot(freq, dbv(spec0[:N/2 + 1]), 'c', linewidth=1, label='$S$')
    plt.hold(True)
    plotSpectrum(spec0, fin, '--r', linewidth=2, label='$\\mathrm{plotSpectrum}(S)$')
    plt.plot(freq, dbp(spec_smoothed[:N/2 + 1]), 'b', linewidth=2, label='$\\mathrm{circ\\_smooth}(S)$')
    plt.plot(freq, dbp(Snn*NBW), 'm', linewidth=1.5, label='$\mathrm{from\\ NTF}$')
    plt.text(0.5, -3, 'NBW = %.1e ' % NBW, horizontalalignment='right',
             verticalalignment='top')
    figureMagic((0, 0.5), None, None, (-140, 0), 20, None)
    ax = plt.gca()
    ax.set_xscale('linear')
    plt.ylabel('Spectrum [dB]')
    plt.xlabel('Normalized frequency ($f_s \\rightarrow 1$)')
    plt.legend(loc=4)
