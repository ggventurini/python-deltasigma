# -*- coding: utf-8 -*-
# _PlotExampleSpectrum.py
# Miscellaneous functions and stdlib wrappers for MATLAB functions
# that do not find a direct replacement in numpy/scipy.
# Copyright 2013 Giuseppe Venturini
# This file is part of python-deltasigma.
#
# python-deltasigma is a 1:1 Python replacement of Richard Schreier's
# MATLAB delta sigma toolbox (aka "delsigma"),  upon which it is heavily based.
# The delta sigma toolbox is (c) 2009,  Richard Schreier.
#
# python-deltasigma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# LICENSE file for the licensing terms.

"""This module provides the PlotExampleSpectrum function.
"""

from __future__ import division
import numpy as np
from numpy.fft import fft, fftshift
import pylab as plt

from ._calculateSNR import calculateSNR
from ._circ_smooth import circ_smooth
from ._dbv import dbv
from ._dbp import dbp
from ._ds_f1f2 import ds_f1f2
from ._ds_hann import ds_hann
from ._evalTF import evalTF
from ._figureMagic import figureMagic
from ._simulateDSM import simulateDSM
from ._undbv import undbv

def PlotExampleSpectrum(ntf, M=1, osr=64, f0=0, quadrature=False):
    """Plot a spectrum suitable to exemplify the NTF performance.

    **Parameters:**

    ntf : scipy 'lti' object, tuple or ndarray
           The first argument may be one of the various supported
           representations for a (SISO) transfer function or an
           ABCD matrix. See :func:`evalTF` for a more detailed
           discussion.

    f0 : float, optional
         The center frequency. Normalized. Defaults to 0.

    M : int, optional
        M is defined as:

        .. math::

            M = n_{lev} - 1

        The number of quantizer levels (:math:`n_{lev}`) defaults to 2 and
        ``M`` defaults to 1.

    quadrature : boolean, optional
                 Whether the delta sigma modulator is a quadrature
                 modulator or not. Defaults to ``False``.

    .. note::

        Quadrature modulators require the following currently unimplemented
        functions: :func:`simulateQDSM` (expected v. 0.2).
        Setting ``quadrature`` to ``True`` results in a
        ``NotImplementedError`` being raised.

    .. plot::

        import pylab as plt
        from deltasigma import synthesizeNTF, PlotExampleSpectrum
        order = 3
        osr = 32
        f0 = 0.
        Hinf = 1.5
        ntf = synthesizeNTF(order, osr, 0, Hinf, f0)
        plt.figure(figsize=(12, 5))
        PlotExampleSpectrum(ntf, M=1, osr=osr, f0=f0)

    """
    f1, f2 = ds_f1f2(osr, f0, quadrature)
    delta = 2
    Amp = undbv(-3) # Test tone amplitude, relative to full-scale.
    # f below is the test tone frequency offset from f0, relative to bw.
    # (It will be adjusted to be an fft bin)
    f = 0.3
    N = 2**12
    f1_bin = np.round(f1*N)
    f2_bin = np.round(f2*N)
    fin = round(((1 - f)/2*f1 + (f + 1)/2*f2) * N)
    if not quadrature:
        t = np.arange(0, N).reshape((1, -1))
        u = Amp*M*np.cos((2*np.pi/N)*fin*t)
        v, _, xmax, y = simulateDSM(u, ntf, M+1)
    else:
        raise NotImplementedError("The required simulateQDSM function " + \
                                  "is not available yet.")
        t = np.arange(0, N).reshape((1, -1))
        u = Amp*M*np.exp((2j*np.pi/N)*fin*t)
        v, xn, xmax, y = simulateQDSM(u, ntf, M + 1)
    window = ds_hann(N)
    NBW = 1.5/N
    spec0 = fft(v * window)/(M*N/4)
    if not quadrature:
        freq = np.linspace(0, 0.5, N/2 + 1)
        plt.plot(freq, dbv(spec0[:N/2 + 1]), 'c', linewidth=1)
        plt.hold(True)
        spec_smoothed = circ_smooth(np.abs(spec0)**2., 16)
        plt.plot(freq, dbp(spec_smoothed[:N/2 + 1]), 'b', linewidth=3)
        Snn = np.abs(evalTF(ntf, np.exp(2j*np.pi*freq)))**2 * 2/12*(delta/M)**2
        plt.plot(freq, dbp(Snn*NBW), 'm', linewidth=1)
        snr = calculateSNR(spec0[f1_bin:f2_bin + 1], fin - f1_bin)
        msg = 'SQNR  =  %.1fdB\n @ A = %.1fdBFS & osr = %.0f\n' % \
              (snr, dbv(spec0[fin]), osr)
        if f0 < 0.25:
            plt.text(f0 + 1 / osr, - 15, msg, horizontalalignment='left',
                     verticalalignment='center')
        else:
            plt.text(f0 - 1 / osr, - 15, msg, horizontalalignment='right',
                     verticalalignment='center')
        plt.text(0.5, - 135, 'NBW = %.1e ' % NBW, horizontalalignment='right',
                 verticalalignment='bottom')
        figureMagic((0, 0.5), 1./16, None, (-140, 0), 10, None)
    else:
        spec0 = fftshift(spec0 / 2)
        freq = np.linspace(-0.5, 0.5, N + 1)
        freq = freq[:-1]
        plt.plot(freq, dbv(spec0), 'c', linewidth=1)
        plt.hold('on')
        spec_smoothed = circ_smooth(abs(spec0) ** 2, 16)
        plt.plot(freq, dbp(spec_smoothed), 'b', linewidth=3)
        Snn = abs(evalTF(ntf, np.exp(2j * np.pi * freq))) ** 2 * 2 / 12 * (delta / M) ** 2
        plt.plot(freq, dbp(Snn * NBW), 'm', linewidth=1)
        snr = calculateSNR(spec0[N/2 + np.arange(f1_bin, f2_bin + 1)], fin - f1_bin)
        msg = 'SQNR  =  %.1fdB\\n @ A = %.1fdBFS & osr = %.0f\\n' % \
              (snr, dbv(spec0[N/2 + fin]), osr)
        if f0 >=  0:
            plt.text(f0 - 0.05, - 15, msg, horizontalalignment='right',
                     verticalalignment='bottom')
        else:
            plt.text(f0 + 0.05, - 15, msg, horizontalalignment='left',
                     verticalalignment='bottom')
        plt.text(-0.5, -135, ' NBW = %.1e' % NBW, horizontalalignment='left',
                 verticalalignment='bottom')
        figureMagic((-0.5, 0.5), 0.125, None, (-140, 0), 10, None)
    plt.xlabel('frequency')
    return

