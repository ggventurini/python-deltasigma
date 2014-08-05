# -*- coding: utf-8 -*-
# test_evalTFP.py
# This module provides the tests for the evalTFP function.
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

"""This module provides the test class for the evalTFP() function.
"""

import unittest
import numpy as np
import deltasigma as ds

class TestEvalTFP(unittest.TestCase):
    """Test class for evalTFP()"""
    def setUp(self):
        pass

    def test_evalTFP_1(self):
        """Test function for evalTFP() 1/2"""
        #          (z-0.3)
        # H1 = ---------------
        #      (z-0.5) (z-0.9)
        H1 = (np.array([.3]), np.array([.5, .9]), 1)
        #             1
        # H2 = ---------------
        #      (z-0.3) (z-0.9)
        H2 = (np.array([]), np.array([.3, .9]), 1)
        a = ds.evalTFP(H1, H2, .2)
        at = np.array([0.5611 + 0.1483j])
        self.assertTrue(np.allclose(np.array([a]), at, atol=1e-4, rtol=1e-4))
        self.assertTrue(np.isscalar(a))

    def test_evalTFP_2(self):
        """Test function for evalTFP() 2/2"""
        #          (z-0.301)
        # H1 = ---------------
        #      (z-0.5) (z-0.9)
        H1 = (np.array([.301]), np.array([.5, .9]), 1)
        #             1
        # H2 = ---------------
        #      (z-0.3) (z-0.9)
        H2 = (np.array([]), np.array([.3, .9]), 1)
        a = ds.evalTFP(H1, H2, np.array([.2, .23, .4]))
        at = np.array([0.5610 + 0.1488j, 0.4466 + 0.0218j, 0.0632 - 0.1504j])
        self.assertTrue(np.allclose(a, at, atol=1e-4, rtol=1e-4))
        self.assertTrue((3,) == a.shape)

