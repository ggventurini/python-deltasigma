
# -*- coding: utf-8 -*-
# test_rms.py
# This module provides the tests for the rms function.
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

"""This module provides the test class for the rms() function.
"""

import unittest
import numpy as np

from deltasigma import rms

class TestRms(unittest.TestCase):
    """Test class for rms()"""

    def setUp(self):
        self.tv = np.arange(100)
        self.res1 = np.sqrt(np.sum(self.tv**2.)/float(self.tv.shape[0]))
        self.res2 = np.sqrt((np.sum((self.tv - self.tv.mean())**2.)) \
                    /self.tv.shape[0])

    def test_rms_1(self):
        """Test function for rms() 1/2"""
        self.assertTrue(np.allclose(rms(self.tv), self.res1, rtol=1e-05,
                                    atol=1e-08))

    def test_rms_2(self):
        """Test function for rms() 2/2"""
        self.assertTrue(np.allclose(rms(self.tv, no_dc=True), self.res2,
                                    rtol=1e-05, atol=1e-08))

