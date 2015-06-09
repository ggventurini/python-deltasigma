# -*- coding: utf-8 -*-
# test_simulateQSNR.py
# This module provides the tests for the simulateQSNR function.
# Copyright 2015 Giuseppe Venturini
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

"""This module provides the test class for the simulateSNR() function.
"""

from __future__ import division, print_function

import unittest
import numpy as np
import deltasigma as ds

class TestSimulateQSNR(unittest.TestCase):
    """Test class for simulateQSNR()"""

    def setUp(self):
        self.SNR_ref = np.array([23.0421, 32.1100, 43.3758, 53.1791,
                                 65.5504, 70.5023, 73.4608, 76.2416, 77.8770,
                                 78.2733, 79.3729, 79.5728, 80.8729, 82.7461,
                                 83.0723, 84.8488, 84.3327])
        self.AMP_ref = np.array([-70, -60, -50, -40, -30, -20, -15, -10, -9, -8,
                                 -7, -6, -5, -4, -3, -2, -1, 0])


    def test_simulateQSNR(self):
        """Test function for simulateSNR() 4/4"""
        order = 4
        osr = 32
        M = 8
        NG = -50
        ING = -10
        f0 = 1./ 16
        quadrature = 1
        form = 'PFB'
        nlev = M + 1
        z0 = np.exp(1j*2*np.pi*f0)
        bw = 1./ osr
        delta = 2
        FullScale = M
        ntf0 = ds.synthesizeQNTF(order, osr, f0, NG, ING)
        ABCD = ds.realizeQNTF(ntf0, form, True)
        a, b = ds.simulateQSNR(ABCD, osr, None, f0, nlev);
        assert np.allclose(a[6:], self.SNR_ref, atol=6, rtol=1e-3)
        assert np.allclose(b[5:], self.AMP_ref, atol=1)
