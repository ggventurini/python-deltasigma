# -*- coding: utf-8 -*-
# test_impL1.py
# Test module for the impL1() function.
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

"""Test module for the impL1() function.
"""

import unittest
import numpy as np
import deltasigma as ds

from scipy.signal import lti

class TestTestTemplate(unittest.TestCase):
    """Test class for the impL1 function"""
    def setUp(self):
        self.r2 = np.array([0., 0.4, -0.16, 0.064, -0.0256, 0.01024,
                            -0.004096, 0.0016384, -0.00065536,
                            0.000262144, -0.0001048576])

    def test_impL1_zpk(self):
        """Test function for impL1() 1/4"""
        sys1 = (np.array([-.4]), np.array([0, 0]), 1) # zpk
        r1 = ds.impL1(sys1, n=10)
        self.assertTrue(np.allclose(r1, self.r2, atol=1e-8, rtol=1e-4))

    def test_impL1_num_den(self):
        """Test function for impL1() 2/4"""
        sys1 = (np.array([-.4]), np.array([0, 0]), 1) 
        num = np.poly(sys1[0])
        den = np.poly(sys1[1])
        num[0] *= sys1[2]
        tf = (num, den)
        r3 = ds.impL1(tf, n=10)
        self.assertTrue(np.allclose(self.r2, r3, atol=1e-8, rtol=1e-4))

    def test_impL1_lti(self):
        """Test function for impL1() 3/4"""
        sys1 = (np.array([-.4]), np.array([0, 0]), 1)
        tf = lti(*sys1)
        r4 = ds.impL1(tf, n=10)
        self.assertTrue(np.allclose(self.r2, r4, atol=1e-8, rtol=1e-4))

    def test_impL1_cplx(self):
        """Test function for impL1() 4/4"""
        res = np.array([0.00000000000000 + 0.00000000000000j,
                        0.333000000000000 - 0.333000000000000j,
                        -0.110889000000000 + 0.110889000000000j,
                        0.0369260370000000 - 0.0369260370000000j,
                        -0.0122963703210000 + 0.0122963703210000j,
                        0.00409469131689300 - 0.00409469131689300j,
                        -0.00136353220852537 + 0.00136353220852537j,
                        0.000454056225438948 - 0.000454056225438948j,
                        -0.000151200723071170 + 0.000151200723071170j,
                        5.03498407826996e-05 - 5.03498407826996e-05j,
                        -1.67664969806390e-05 + 1.67664969806390e-05j])

        a = [[-.3333], [-.3333j], 1]
        y = ds.impL1(a)
        self.assertTrue(np.allclose(y, res, atol=1e-3, rtol=1e-4)) 
