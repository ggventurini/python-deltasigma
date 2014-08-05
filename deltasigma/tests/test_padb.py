# -*- coding: utf-8 -*-
# test_padb.py
# This module provides the tests for the padb function.
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

"""This module provides the test class for the padb() function.
"""

import unittest
import numpy as np
from deltasigma import padb

class TestPadb(unittest.TestCase):
    """Test class for padb()"""
    def setUp(self):
        self.tv1 = np.eye(15)
        self.tv2 = np.arange(10)
        self.tv3 = np.array([])

    def test_padb_1(self):
        """Test function for padb() 1/3"""
        # regular matrix
        tr = padb(self.tv1, n=25, val=2)
        res = np.concatenate((self.tv1, 2.*np.ones((10, 15))), axis=0)
        self.assertTrue(np.allclose(tr, res, atol=1e-8, rtol=1e-5))

    def test_padb_2(self):
        """Test function for padb() 2/3"""
        # 1-d array
        tr = padb(self.tv2, n=25, val=1.5)
        res = np.vstack((self.tv2.reshape((-1, 1)), 1.5*np.ones((15, 1))))
        self.assertTrue(np.allclose(tr, res, atol=1e-8, rtol=1e-5))

    def test_padb_3(self):
        """Test function for padb() 3/3"""
        # empty matrix array
        tr = padb(self.tv3, n=25, val=1.5)
        res = np.vstack((self.tv3.reshape((-1, 1)), 1.5*np.ones((25, 1))))
        self.assertTrue(np.allclose(tr, res, atol=1e-8, rtol=1e-5))

