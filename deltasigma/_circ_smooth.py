# -*- coding: utf-8 -*-
# _circ_smooth.py
# Module providing the circ_smooth function
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

"""Module providing the circ_smooth() function
"""

from __future__ import division
import numpy as np
from ._ds_hann import ds_hann
from ._utils import circshift

def circ_smooth(x, n=16):
    """Smoothing of the PSD ``x`` for linear x-axis plotting.

    **Parameters:**

    x : 1D ndarray
        The PSD to be smoothed, or equivalently

    .. math::

        x = \\left|\\mathrm{FFT}\\left(u(t)\\right)(f)\\right|^2

    n : int, even, optional
        The length of the Hann window used in the smoothing algorithm.

    **Returns:**

    y : 1D ndarray
        The smoothed PSD

    .. seealso::

        :func:`logsmooth`, smoothing algorithm suitable for logarithmic x-axis
        plotting.

    For a comparison of :func:`circ_smooth` and :func:`logsmooth` (accessed
    through the helper function :func:`plotSpectrum`) see the following plot.

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
        # plotting
        plt.figure(figsize=(12, 7))
        plt.subplot(211)
        plt.plot(freq, dbv(spec0[:N/2 + 1]), 'c', linewidth=1, label='$S$')
        plt.hold(True)
        spec_smoothed = circ_smooth(np.abs(spec0)**2., 16)
        plt.plot(freq, dbp(spec_smoothed[:N/2 + 1]), 'b--', linewidth=2, label='$\\\\mathrm{circ\\\\_smooth}(S)$')
        plotSpectrum(spec0, fin, 'r', linewidth=2, label='$\\\\mathrm{plotSpectrum}(S)$')
        Snn = np.abs(evalTF(ntf, np.exp(2j*np.pi*freq)))**2 * 2/12*(delta)**2
        plt.plot(freq, dbp(Snn*NBW), 'm', linewidth=1.5, label='$\\\\mathrm{from\\\\ NTF}$')
        snr = calculateSNR(spec0[f1_bin:f2_bin + 1], fin - f1_bin)
        msg = 'SQNR  =  %.1fdB\\n @ A = %.1fdBFS & osr = %.0f\\n' % \\
              (snr, dbv(spec0[fin]), osr)
        plt.text(f0 + .45/osr, -20, msg, horizontalalignment='left',
                 verticalalignment='center')
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
        plotSpectrum(spec0, fin, '--r', linewidth=2, label='$\\\\mathrm{plotSpectrum}(S)$')
        plt.plot(freq, dbp(spec_smoothed[:N/2 + 1]), 'b', linewidth=2, label='$\\\\mathrm{circ\\\\_smooth}(S)$')
        plt.plot(freq, dbp(Snn*NBW), 'm', linewidth=1.5, label='$\\\\mathrm{from\\\\ NTF}$')
        msg = 'SQNR  =  %.1fdB\\n @ A = %.1fdBFS & osr = %.0f\\n' % \\
              (snr, dbv(spec0[fin]), osr)
        plt.text(f0 + 1./osr, -20, msg, horizontalalignment='left',
                 verticalalignment='center')
        plt.text(0.5, -3, 'NBW = %.1e ' % NBW, horizontalalignment='right',
                 verticalalignment='top')
        figureMagic((0, 0.5), None, None, (-140, 0), 20, None)
        ax = plt.gca()
        ax.set_xscale('linear')
        plt.ylabel('Spectrum [dB]')
        plt.xlabel('Normalized frequency ($f_s \\\\rightarrow 1$)')
        plt.legend(loc=4)

    """
    assert len(x.shape) == 1 or 1 in x.shape
    assert n % 2 == 0
    nx = max(x.shape)
    w = ds_hann(n)/(n/2.)
    xw = np.convolve(x, w)
    yp = np.hstack((xw[n - 1:nx], xw[:n - 1] + xw[nx:]))
    y = circshift(yp, [int(n/2. - 1)])
    return y

