# -*- coding: utf-8 -*-
# test_bunquantize.py
# This module provides the tests for the bunquantize function.
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

"""This module provides the test class for the bunquantize() function.
"""

import unittest
import numpy as np
import deltasigma as ds

class TestBUnQuantize(unittest.TestCase):
    """Test class for bunquantize()"""
    def setUp(self):
        self.x = np.linspace(-10, 10, 101)

    def test_bunquantize(self):
        """Test function for bunquantize() 1/3"""
        yr = ds.bquantize(self.x)
        yv = []
        y = []
        for yi in yr:
            y += [yi.csd]
            yv += [yi.val]
        yv = np.asarray(yv)
        xres = ds.bunquantize(y)
        self.assertTrue(np.allclose(xres, yv, atol=1e-8, rtol=1e-5))

