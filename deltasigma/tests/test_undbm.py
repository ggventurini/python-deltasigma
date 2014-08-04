# -*- coding: utf-8 -*-
# test_undbm.py
# This module provides the tests for the undbm function.
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

"""This module provides the test class for the undbm() function.
"""

import unittest
import numpy as np

from deltasigma import undbm

class TestUndbm(unittest.TestCase):
    """Test class for undbm()"""

    def setUp(self):
        pass

    def test_undbm_1(self):
        """Test function for undbm() 1/3"""
        self.assertTrue(np.allclose([undbm(53.015)], [100.054125892], rtol=1e-05,
                                    atol=1e-08))

    def test_undbm_2(self):
        """Test function for undbm() 2/3"""
        self.assertTrue(np.allclose([undbm(3, 100)], [0.44668359215], rtol=1e-05,
                                    atol=1e-08))

    def test_undbm_3(self):
        """Test function for undbm() 3/3"""
        self.assertTrue(np.isscalar(undbm(3, 100)))

