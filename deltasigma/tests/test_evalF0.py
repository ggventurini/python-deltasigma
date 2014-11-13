# -*- coding: utf-8 -*-
# test_evalF0.py
# This module provides the tests for the evalF0 function.
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

"""This module provides the test class for the evalF0() function.
"""

import unittest
import numpy as np
from math import pi

from deltasigma import evalF0

class TestEvalF0(unittest.TestCase):
    """Test class for evalF0()"""

    def setUp(self):
        self.z = np.exp(1j*2*pi*np.arange(11.0)/10.0)

    def test_evalF0_1(self):
        """Test function for evalF0() 1/2"""
        r = evalF0(1, self.z, 1)
        ref = np.array((1.5, 1.309016994374948, 0.809016994374947,
                        0.190983005625053, -0.309016994374947, -0.5,
                        -0.309016994374947, 0.190983005625052,
                        0.809016994374947, 1.309016994374947, 1.5))
        self.assertTrue(np.allclose(r, ref))

    def test_evalF0_2(self):
        """Test function for evalF0() 2/3"""
        r = evalF0(1, self.z, 0.5)
        ref = np.array((2.5, 2.118033988749895, 1.118033988749895,
                        -0.118033988749895, -1.118033988749895, -1.5,
                        -1.118033988749895, -0.118033988749895,
                        1.118033988749895, 2.118033988749895, 2.5))
        self.assertTrue(np.allclose(r, ref))

    def test_evalF0_3(self):
        """Test function for evalF0() 3/3"""
        r = evalF0(0.5, self.z, 0.5)
        ref = np.array((1.5, 1.309016994374948, 0.809016994374947,
                        0.190983005625053, -0.309016994374947, -0.5,
                        -0.309016994374947,  0.190983005625052,
                        0.809016994374947, 1.309016994374947, 1.5))
        self.assertTrue(np.allclose(r, ref))

