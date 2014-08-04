# -*- coding: utf-8 -*-
# test_undbv.py
# This module provides the tests for the undbv function.
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

"""This module provides the test class for the undbv() function.
"""

import unittest
import numpy as np

from deltasigma import undbv

class TestUndbv(unittest.TestCase):
    """Test class for undbv()"""

    def setUp(self):
        pass

    def test_undbv_1(self):
        """Test function for undbv() 1/2"""
        self.assertTrue(np.allclose([undbv(53.05)], [449.26232467], rtol=1e-05,
                                    atol=1e-08))

    def test_undbv_2(self):
        """Test function for undbv() 2/2"""
        self.assertTrue(np.allclose([undbv(3)], [1.41253754462], rtol=1e-05,
                                    atol=1e-08))
        


