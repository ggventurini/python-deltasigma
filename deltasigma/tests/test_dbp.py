# -*- coding: utf-8 -*-
# test_dbp.py
# This module provides the tests for the dbp function.
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

"""This module provides the test class for the dbp() function.
"""

import unittest
import numpy as np

from deltasigma import dbp

class TestDbp(unittest.TestCase):
    """Test class for dbp()"""
    def setUp(self):
        self.tv1 = np.array([2])
        self.r1 = np.array([3.01029996])
        self.tv2 = 2
        self.r2 = 3.01029996
        self.tv3 = 2, 2
        self.r3 = 3.01029996, 3.01029996

    def test_dbp_1(self):
        """Test function for dbp() 1/3"""
        res = dbp(self.tv1)
        self.assertTrue(np.allclose(self.r1, res, atol=1e-8, rtol=1e-5))

    def test_dbp_2(self):
        """Test function for dbp() 2/3"""
        res = dbp(self.tv2)
        self.assertTrue(np.allclose(self.r2, res, atol=1e-8, rtol=1e-5))
        self.assertTrue(np.isscalar(res))  # check for type coherence

    def test_dbp_3(self):
        """Test function for dbp() 3/3"""
        res = dbp(self.tv3)
        self.assertTrue(np.allclose(self.r3, res, atol=1e-8, rtol=1e-5))

