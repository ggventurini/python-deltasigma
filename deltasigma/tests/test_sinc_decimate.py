# -*- coding: utf-8 -*-
# test_sinc_decimate.py
# This module provides the tests for the sinc_decimate function.
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

"""This module provides the test class for the sinc_decimate() function.
"""

import unittest
import numpy as np
import deltasigma as ds

class TestSincDecimate(unittest.TestCase):
    """Test class for sinc_decimate()"""

    def test_sinc_decimate_1(self):
        """Test function for sinc_decimate() 1/3"""
        x = [1]*10
        self.assertTrue(np.allclose(ds.sinc_decimate(x, 1, 10), [1.]))

    def test_sinc_decimate_2(self):
        """Test function for sinc_decimate() 2/3"""
        x = [1]*10
        self.assertTrue(np.allclose(ds.sinc_decimate(x, 2, 5), [.6, 1.]))

    def test_sinc_decimate_3(self):
        """Test function for sinc_decimate() 3/3"""
        x = [1]*10
        x = np.cumsum(np.asarray(x))
        self.assertTrue(np.allclose(ds.sinc_decimate(x, 1, 5), [3., 8.]))
        self.assertTrue(np.allclose(ds.sinc_decimate(x, 2, 5), [1.4, 6.]))
        self.assertTrue(np.allclose(ds.sinc_decimate(x, 3, 3), [0.55555556, 3., 6.]))

