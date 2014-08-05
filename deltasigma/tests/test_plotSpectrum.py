# -*- coding: utf-8 -*-
# test_plotSpectrum.py
# This module provides the tests for the plotSpectrum function.
# Copyright 2014 Giuseppe Venturini
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

"""This module provides the test class for the plotSpectrum() function.
"""

import unittest
import numpy as np
import deltasigma as ds
import pylab as plt

from numpy.fft import fft

class TestPlotSpectrum(unittest.TestCase):
    """Test class for plotSpectrum()"""

    def setUp(self):
        pass

    def test_plotSpectrum(self):
        """Test function for plotSpectrum()"""
        f0 = 0
        osr = 32
        quadrature = False
        Hinf = 1.5
        order = 3
        ntf = ds.synthesizeNTF(order, osr, 0, Hinf, f0)
        f1, f2 = ds.ds_f1f2(osr, f0, quadrature)
        delta = 2
        Amp = ds.undbv(-3)
        f = 0.3
        N = 2**12
        f1_bin = np.round(f1*N)
        f2_bin = np.round(f2*N)
        fin = np.round(((1 - f)/2*f1 + (f + 1)/2*f2) * N)
        t = np.arange(0, N)
        u = Amp*np.cos((2*np.pi/N)*fin*t)
        v, xn, xmax, y = ds.simulateDSM(u, ntf, 2)
        window = ds.ds_hann(N)
        NBW = 1.5/N
        spec0 = fft(v * window)/(N/4)
        freq = np.linspace(0, 0.5, N/2 + 1)
        # plotting
        plt.subplot(211)
        plt.plot(freq, ds.dbv(spec0[:N/2 + 1]), 'c', linewidth=1, label='$S$')
        plt.hold(True)
        spec_smoothed = ds.circ_smooth(np.abs(spec0)**2., 16)
        plt.plot(freq, ds.dbp(spec_smoothed[:N/2 + 1]), 'b--', linewidth=2,
                 label='$\\mathrm{circ\\_smooth}(S)$')
        ds.plotSpectrum(spec0, fin, 'r', linewidth=2,
                        label='$\\mathrm{plotSpectrum}(S)$')
        Snn = np.abs(ds.evalTF(ntf, np.exp(2j*np.pi*freq)))**2 * 2/12*(delta)**2
        plt.plot(freq, ds.dbp(Snn*NBW), 'm', linewidth=1.5,
                 label='$\mathrm{from\\ NTF}$')
        plt.text(0.5, -3, 'NBW = %.1e ' % NBW, horizontalalignment='right',
                 verticalalignment='top')
        ds.figureMagic((0, 0.5), None, None, (-140, 0), 20, None)
        plt.ylabel('Spectrum [dB]')
        ax = plt.gca()
        ax.set_title('Smoothing and plotting for LOG and LIN axes')
        plt.legend(loc=4)
        plt.subplot(212)
        plt.plot(freq, ds.dbv(spec0[:N/2 + 1]), 'c', linewidth=1, label='$S$')
        plt.hold(True)
        ds.plotSpectrum(spec0, fin, '--r', linewidth=2,
                        label='$\\mathrm{plotSpectrum}(S)$')
        plt.plot(freq, ds.dbp(spec_smoothed[:N/2 + 1]), 'b', linewidth=2,
                 label='$\\mathrm{circ\\_smooth}(S)$')
        plt.plot(freq, ds.dbp(Snn*NBW), 'm', linewidth=1.5,
                 label='$\mathrm{from\\ NTF}$')
        plt.text(0.5, -3, 'NBW = %.1e ' % NBW, horizontalalignment='right',
                 verticalalignment='top')
        ds.figureMagic((0, 0.5), None, None, (-140, 0), 20, None)
        ax = plt.gca()
        ax.set_xscale('linear')
        plt.ylabel('Spectrum [dB]')
        plt.xlabel('Normalized frequency ($f_s \\rightarrow 1$)')
        plt.legend(loc=4)

