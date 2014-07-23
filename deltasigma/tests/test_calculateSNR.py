# -*- coding: utf-8 -*-
# test_calculateSNR.py
# This module provides the tests for the calculateSNR() function.
# Copyright 2014 Giuseppe Venturini & Shayne Hodge
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

import unittest
import numpy as np
import deltasigma as ds
from numpy.fft import fft
from deltasigma._ds_hann import ds_hann


class TestCalculateSNR(unittest.TestCase):
    """Test function for calculateSNR()"""
    def setUp(self):
        N = 2**12
        t = np.arange(N)
        (f1, f2) = (1.0/8, 1.0/302)
        A = np.cos(2*np.pi*f1*t)
        B = 0.01*np.cos(2*np.pi*f2*t)
        y = A + B
        window = ds_hann(N)
        self.hwfft = fft(window*y)
        self.N = N
        self.f1 = f1

    def test_snr_is_40(self):
        """ Test that a particular SNR is within roundings errors of
        40 (dB?) """
        N = self.N
        snr = ds.calculateSNR(self.hwfft[:N/2], int(N*self.f1))
        # Consider replacing with assertAlmostEqual
        self.assertTrue(np.allclose(snr, 40, atol=1e-8, rtol=1e-8))

    def test_snr_is_inf(self):
        """ Test that a paricular SNR is infinite. """
        N = self.N
        hwfft = np.zeros((N/2, ))
        hwfft[512] = 1.0  # specially crafted to have Inf snr
        snr = ds.calculateSNR(hwfft[:N/2], 512)
        self.assertEqual(snr, np.Inf)
