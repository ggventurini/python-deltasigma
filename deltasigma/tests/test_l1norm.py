# -*- coding: utf-8 -*-
# test_l1norm.py
# This module provides the tests for the l1norm function.
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

"""This module provides the test class for the l1norm() function.
"""

import unittest
import numpy as np
import deltasigma as ds

from scipy.signal import lti

class Testl1norm(unittest.TestCase):
    """Test class for l1norm()"""

    def setUp(self):
        zeros = np.array(())
        poles = np.array((.5,))
        k = 1.
        self.zpk_tuple = zeros, poles, k
        splti = lti(zeros, poles, k)
        self.num_den_tuple = (splti.num, splti.den)
        self.ABCD_tuple = (splti.A, splti.B, splti.C, splti.D)
        self.splti = splti

    def test_l1norm_1(self):
        """Test function for l1norm() 1/4"""
        self.assertTrue(np.allclose(ds.l1norm(self.num_den_tuple), 2.,
                        rtol=1e-5, atol=1e-8))

    def test_l1norm_2(self):
        """Test function for l1norm() 2/4"""
        self.assertTrue(np.allclose(ds.l1norm(self.zpk_tuple), 2.,
                        rtol=1e-5, atol=1e-8))

    def test_l1norm_3(self):
        """Test function for l1norm() 3/4"""
        self.assertTrue(np.allclose(ds.l1norm(self.ABCD_tuple), 2.,
                        rtol=1e-5, atol=1e-8))

    def test_l1norm_4(self):
        """Test function for l1norm() 4/4"""
        self.assertTrue(np.allclose(ds.l1norm(self.splti), 2.,
                        rtol=1e-5, atol=1e-8))

