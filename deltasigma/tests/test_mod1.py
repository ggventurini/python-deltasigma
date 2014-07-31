# -*- coding: utf-8 -*-
# test_mod1.py
# This module provides the tests for the mod1 function.
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

"""This module provides the test class for the mod1() function.
"""

import unittest
import numpy as np
import deltasigma as ds

class Testmod1(unittest.TestCase):
    """Test class for mod1()"""

    def setUp(self):
        self.ABCDmod1 = [[1., 1., -1.], [1., 0., 0.]]

    def test_mod1(self):
        """Test function for mod1()"""
        ABCD, ntf, stf = ds.mod1()
        self.assertTrue(np.allclose(ABCD, self.ABCDmod1, atol=1e-8, rtol=1e-5))

