# -*- coding: utf-8 -*-
# test_mod2.py
# This module provides the tests for the mod2 function.
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

"""This module provides the test class for the mod2() function.
"""

import unittest
import numpy as np
import deltasigma as ds

class TestMod2(unittest.TestCase):
    """Test class for mod2()"""

    def setUp(self):
        self.ABCDmod2 = [[1., 0., 1., -1.],
                         [1., 1., 1., -2.],
                         [0., 1., 0., 0.]]

    def test_mod2(self):
        """Test function for mod2()"""
        ABCD, ntf, stf = ds.mod2()
        self.assertTrue(np.allclose(ABCD, self.ABCDmod2, atol=1e-8,
                        rtol=1e-5))

