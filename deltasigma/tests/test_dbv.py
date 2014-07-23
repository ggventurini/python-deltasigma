# -*- coding: utf-8 -*-
# test_dbv.py
# This module provides the tests for the dbv function.
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

"""This module provides the test class for the dbv() function.
"""

import unittest
import numpy as np

from deltasigma import dbv, undbv

class TestDbv(unittest.TestCase):
    """Test class for dbv()"""
    def setUp(self):
        self.t1 = np.array([3.0])
        self.r1 = 9.5424250943932485

    def test_dbv(self):
        """Test function for dbv()"""
        res = dbv(self.t1)
        t2 = undbv(res)
        self.assertTrue(np.allclose(t2, self.t1, atol=1e-8, rtol=1e-5))
        self.assertTrue(np.allclose(res, self.r1, atol=1e-8, rtol=1e-5))

