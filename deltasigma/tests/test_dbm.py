# -*- coding: utf-8 -*-
# test_dbm.py
# This module provides the tests for the dbm function.
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

"""This module provides the test class for the dbm() function.
"""

import unittest
import numpy as np

from deltasigma import dbm

class TestDbm(unittest.TestCase):
    """Test class for dbm()"""
    def setUp(self):
        self.v1 = np.arange(10) * 1e-3 # test arrays
        self.r1 = [-np.inf, -46.98970004, -40.96910013, -37.44727495, -34.94850022,
                   -33.01029996, -31.42667504, -30.08773924, -28.9279003, -27.90484985]
        self.v2 = 9e-3  # test scalars.
        self.r2 = -27.90484985

    def test_dbm_1(self):
        """Test function for dbm() 1/2"""
        self.assertTrue(np.allclose(dbm(self.v1), self.r1, atol=1e-8, rtol=1e-5))

    def test_dbm_2(self):
        """Test function for dbm() 2/2"""
        self.assertTrue(np.allclose(dbm(self.v2), self.r2, atol=1e-8, rtol=1e-5))
        self.assertTrue(np.isscalar(dbm(self.v2)))

