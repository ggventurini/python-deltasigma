# -*- coding: utf-8 -*-
# test_mapQtoR.py
# This module provides the tests for the mapQtoR function.
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

"""This module provides the test class for the mapQtoR() function.
"""

from __future__ import division

import unittest
import numpy as np
import deltasigma as ds

class TestMapQtoR(unittest.TestCase):
    """Test class for mapQtoR()"""

    def setUp(self):
        A = np.arange(1, 6*7 + 1, dtype=np.int16).reshape((7, 6)).T
        self.A = A + 1j*A
        self.Ares = \
        np.array([[1, -1, 7, -7, 13, -13, 19, -19, 25, -25, 31, -31, 37, -37],
                  [1, 1, 7, 7, 13, 13, 19, 19, 25, 25, 31, 31, 37, 37],
                  [2, -2, 8, -8, 14, -14, 20, -20, 26, -26, 32, -32, 38, -38],
                  [2, 2, 8, 8, 14, 14, 20, 20, 26, 26, 32, 32, 38, 38],
                  [3, -3, 9, -9, 15, -15, 21, -21, 27, -27, 33, -33, 39, -39],
                  [3, 3, 9, 9, 15, 15, 21, 21, 27, 27, 33, 33, 39, 39],
                  [4, -4, 10, -10, 16, -16, 22, -22, 28, -28, 34, -34, 40, -40],
                  [4, 4, 10, 10, 16, 16, 22, 22, 28, 28, 34, 34, 40, 40],
                  [5, -5, 11, -11, 17, -17, 23, -23, 29, -29, 35, -35, 41, -41],
                  [5, 5, 11, 11, 17, 17, 23, 23, 29, 29, 35, 35, 41, 41],
                  [6, -6, 12, -12, 18, -18, 24, -24, 30, -30, 36, -36, 42, -42],
                  [6, 6, 12, 12, 18, 18, 24, 24, 30, 30, 36, 36, 42, 42]],
                 dtype=np.int16)

    def test_mapQtoR(self):
        """Test function for mapQtoR()"""
        At = ds.mapQtoR(self.A)
        self.assertTrue(np.allclose(At, self.Ares))

