# -*- coding: utf-8 -*-
# test_zinc.py
# This module provides the tests for the zinc function.
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

"""This module provides the test class for the zinc() function.
"""

import unittest
import numpy as np
import deltasigma as ds

class TestZinc(unittest.TestCase):
    """Test class for zinc()"""

    def setUp(self):
        self.ref = [1.0000, 0.9985, 0.9941, 0.9867, 0.9765, 0.9635,
                    0.9478, 0.9295, 0.9087, 0.8855, 0.8602, 0.8329,
                    0.8038, 0.7730, 0.7408, 0.7074, 0.6729, 0.6377,
                    0.6019, 0.5658, 0.5295, 0.4933, 0.4574, 0.4221,
                    0.3874, 0.3536, 0.3208, 0.2892, 0.2590, 0.2302,
                    0.2031, 0.1776, 0.1538, 0.1319, 0.1118, 0.0936,
                    0.0772, 0.0626, 0.0499, 0.0389, 0.0295, 0.0217,
                    0.0154, 0.0104, 0.0066, 0.0038, 0.0020, 0.0008,
                    0.0002, 0.0000, 0.0000]
        self.f = np.arange(0, 0.51, 0.01)

    def test_zinc(self):
        """Test function for zinc()"""
        test = ds.zinc(self.f, 2, 3)
        self.assertTrue(np.allclose(test, self.ref, atol=1e-4,
                                    rtol=1e-4))

