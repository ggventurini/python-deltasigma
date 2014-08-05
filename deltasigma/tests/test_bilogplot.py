# -*- coding: utf-8 -*-
# test_bilogplot.py
# Test module for the bilogplot function
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

"""Test module for the bilogplot function
"""

from __future__ import division

import unittest
import numpy as np
import pylab as plt
import deltasigma as ds

class TestBiLogPlot(unittest.TestCase):
    """Test function for bilogplot()"""

    def setUp(self):
        pass

    def test_bilogplot(self):
        """Test function for bilogplot()"""
        f0 = 1./8
        OSR = 64
        order = 8
        N = 8192
        H = ds.synthesizeNTF(order, OSR, 1, 1.5, f0)
        fB = int(np.ceil(N/(2. * OSR)))
        ftest = int(np.round(f0*N + 1./3 * fB))
        u = 0.5*np.sin(2*np.pi*ftest/N*np.arange(N))
        v, xn, xmax, y = ds.simulateDSM(u, H)
        spec = np.fft.fft(v*ds.ds_hann(N))/(N/4)
        X = spec[:N/2 + 1]
        plt.figure()
        # graphical function: we check it doesn't fail
        ds.bilogplot(X, f0*N, ftest, (.03, .3, .3), (-140, 0, 10, 20))
