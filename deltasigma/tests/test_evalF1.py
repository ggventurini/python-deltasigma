# -*- coding: utf-8 -*-
# test_evalF1.py
# This module provides the tests for the evalF1 function.
# Copyright 2014 Giuseppe Venturini and Shayne Hodge
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

"""This module provides the test class for the evalF1() function.
"""

import unittest
import numpy as np

from deltasigma import evalF1

class TestEvalF1(unittest.TestCase):
    """Test class for evalF1()"""

    def setUp(self):
        pass

    def test_evalF1_1(self):
        """Test function for evalF1() 1/2"""
        r = evalF1([0.5, 1, 1.5, 2, 5, 10, 20, 30, 40.7, 50], 2, 23)
        ref = 0.544143311383570
        self.assertTrue(np.allclose((r,), (ref,)))

    def test_evalF1_2(self):
        """Test function for evalF1() 2/2"""
        r = evalF1([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 0.135)
        ref = 0.640058615996223
        self.assertTrue(np.allclose((r,), (ref,)))

