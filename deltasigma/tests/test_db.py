# -*- coding: utf-8 -*-
# test_db.py
# This module provides the tests for the db function.
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

"""This module provides the test class for the db() function.
"""

import unittest
import numpy as np
import deltasigma as ds

from nose.tools import raises

class TestDB(unittest.TestCase):
    """Test class for db()"""

    def setUp(self):
        # arrays
        self.tv1 = np.array([2])
        self.r1 = np.array([3.01029996])
        # scalars
        self.tv2 = 2
        self.r2 = 3.01029996
        # tuples
        self.tv3 = 2, 2
        self.r3 = 3.01029996, 3.01029996
        # db voltage - undbv
        self.tv4 = np.array([3.0])
        pass

    def test_db_1(self):
        """Test function for db() 1/6"""
        res = ds.db(self.tv1, 'power')
        self.assertTrue(np.allclose(self.r1, res, atol=1e-8, rtol=1e-5))

    def test_db_2(self):
        """Test function for db() 2/6"""
        res = ds.db(self.tv2, 'power')
        self.assertTrue(np.allclose(self.r2, res, atol=1e-8, rtol=1e-5))

    def test_db_3(self):
        """Test function for db() 3/6"""
        res = ds.db(self.tv3, 'power')
        self.assertTrue(np.allclose(self.r3, res, atol=1e-8, rtol=1e-5))

    def test_db_4(self):
        """Test function for db() 4/6"""
        res = ds.undbv(ds.db(self.tv4, 'voltage'))
        self.assertTrue(np.allclose(self.tv4, res, atol=1e-8, rtol=1e-5))

    @raises(ValueError)
    def test_db_5(self):
        """Test function for db() 5/6"""
        res = ds.db([1], 'wrong')

    def test_db_6(self):
        """Test function for db() 6/6"""
        res = ds.db(self.tv1, 'power', R=100)
        # the R value should be ignored (warning only)
        self.assertTrue(np.allclose(self.r1, res, atol=1e-8, rtol=1e-5))

