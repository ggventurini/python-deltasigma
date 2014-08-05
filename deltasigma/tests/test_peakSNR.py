# -*- coding: utf-8 -*-
# test_peakSNR.py
# This module provides the tests for the peakSNR function.
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

"""This module provides the test class for the peakSNR() function.
"""

import unittest
import pkg_resources

import numpy as np
import scipy.io

import deltasigma as ds

class TestPeakSNR(unittest.TestCase):
    """Test class for peakSNR()"""

    def setUp(self):
        ds._peakSNR._debug = True
        fname = pkg_resources.resource_filename(__name__, "test_data/test_peak_snr.mat")
        self.snr = scipy.io.loadmat(fname)['snr'].reshape((-1,))
        self.amp = scipy.io.loadmat(fname)['amp'].reshape((-1,))
        self.peak_snr, self.peak_amp = 76.612340603949761, -3.220409771005124

    def test_peakSNR(self):
        """Test function for peakSNR()"""
        ps, pa = ds.peakSNR(self.snr, self.amp)
        self.assertTrue(np.allclose(ps, self.peak_snr, atol=1e-8, rtol=1e-5))
        self.assertTrue(np.allclose(pa, self.peak_amp, atol=1e-8, rtol=1e-5))

