# -*- coding: utf-8 -*-
# test_padt.py
# This module provides the tests for the padt function.
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

"""This module provides the test class for the padt() function.
"""

import unittest
import numpy as np

from deltasigma import padt

class TestPadt(unittest.TestCase):
    """Test class for padt()"""
    def setUp(self):
        self.tv1 = np.eye(15)
        self.tv2 = np.arange(10)
        self.tv3 = np.array([])

    def test_padt_1(self):
        """Test function for padt() 1/3"""
        tr = padt(self.tv1, n=25, val=2)
        res = np.concatenate((2.*np.ones((10, 15)), self.tv1), axis=0)
        self.assertTrue(np.allclose(tr, res, atol=1e-8, rtol=1e-5))

    def test_padt_2(self):
        """Test function for padt() 2/3"""
        # 1-d array
        tr = padt(self.tv2, n=25, val=1.5)
        res = np.vstack((1.5*np.ones((15, 1)), self.tv2.reshape((-1, 1))))
        self.assertTrue(np.allclose(tr, res, atol=1e-8, rtol=1e-5))

    def test_padt_3(self):
        """Test function for padt() 3/3"""
        # empty matrix array
        tr = padt(self.tv3, n=25, val=1.5)
        res = np.vstack((1.5*np.ones((25, 1)), self.tv3.reshape((-1, 1))))
        self.assertTrue(np.allclose(tr, res, atol=1e-8, rtol=1e-5))

