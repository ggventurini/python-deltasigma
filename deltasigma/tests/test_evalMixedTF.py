# -*- coding: utf-8 -*-
# test_evalMixedTF.py
# This module provides the tests for the evalMixedTF function.
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

"""This module provides the test class for the evalMixedTF() function.
"""

import unittest
import numpy as np
import deltasigma as ds

class TestEvalMixedTF(unittest.TestCase):
    """Test class for evalMixedTF()"""
    def setUp(self):
        pass

    def test_evalMixedTF_1(self):
        """Test function for evalMixedTF() 1/3"""
        #          (z-0.3)
        # H1 = ---------------
        #      (z-0.5) (z-0.9)
        H1 = (np.array([.3]), np.array([.5, .9]), 1)
        #             1
        # H2 = ---------------
        #      (z-0.3) (z-0.9)
        H2 = (np.array([]), np.array([.3, .9]), 1)
        H = {'Hs':(H1,), 'Hz':(H2,)}
        a = ds.evalMixedTF(H, .2)
        at = np.array([0.5611 + 0.1483j])
        self.assertTrue(np.allclose(np.array([a]), at, atol=1e-4, rtol=1e-4))
        self.assertTrue(np.isscalar(a))

    def test_evalMixedTF_2(self):
        """Test function for evalMixedTF() 2/3"""
        #          (z-0.301)
        # H1 = ---------------
        #      (z-0.5) (z-0.9)
        H1 = (np.array([.301]), np.array([.5, .9]), 1)
        #             1
        # H2 = ---------------
        #      (z-0.3) (z-0.9)
        H2 = (np.array([]), np.array([.3, .9]), 1)
        H = {'Hs':(H1,), 'Hz':(H2,)}
        a = ds.evalMixedTF(H, np.array([.2, .23, .4]))
        at = np.array([0.5610 + 0.1488j, 0.4466 + 0.0218j, 0.0632 - 0.1504j])
        self.assertTrue(np.allclose(a, at, atol=1e-4, rtol=1e-4))
        self.assertTrue((3,) == a.shape)

    def test_evalMixedTF_3(self):
        """Test function for evalMixedTF() 3/3"""
        H11 = (np.array([.3]), np.array([.5, .9]), 1)
        H12 = (np.array([.301]), np.array([.5, .9]), 1)
        H21 = (np.array([]), np.array([.3, .9]), 1)
        H22 = (np.array([]), np.array([.3, .9]), 1)
        H = {'Hs':(H11, H12), 'Hz':(H21, H22)}
        a = ds.evalMixedTF(H, .2)
        at = np.array([0.5611 + 0.1483j])
        self.assertTrue(np.allclose(np.array([a]), 2*at, atol=1e-3, rtol=1e-4))
        self.assertTrue(np.isscalar(a))

