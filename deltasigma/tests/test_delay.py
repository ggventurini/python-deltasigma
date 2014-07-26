# -*- coding: utf-8 -*-
# test_delay.py
# This module provides the tests for the delay function.
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

"""This module provides the test class for the delay() function.
"""

import unittest
import numpy as np

from deltasigma import delay

class TestDelay(unittest.TestCase):
    """Test class for delay()"""

    def setUp(self):
        pass

    def test_delay_1(self):
        """Test function for delay() 1/3"""
        v1 = np.arange(4)
        r1 = np.array((0, 0, 0, 0))
        self.assertTrue(np.allclose(delay(v1, 5), r1))

    def test_delay_2(self):
        """Test function for delay() 2/3"""
        v2 = np.ones((10, 1))
        r2 = np.zeros((10, 1))
        r2[5:, 0] = 1
        self.assertTrue(np.allclose(delay(v2, 5), r2))

    def test_delay_3(self):
        """Test function for delay() 3/3"""
        v3 = np.ones((1, 10))
        r3 = np.zeros((1, 10))
        r3[0, 5:] = 1
        self.assertTrue(np.allclose(delay(v3, 5), r3))

