# -*- coding: utf-8 -*-
# test_cancelPZ.py
# This module provides the tests for the cancelPZ() function.
# Copyright 2014 Giuseppe Venturini & Shayne Hodge
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

import unittest
import numpy as np
import deltasigma as ds


class TestCancelPZ(unittest.TestCase):
    """Test function for calculateTF()"""
    def setUp(self):
        zt = np.array((1, 2, 3))
        pt = np.array((1 + 2e-5, 2 + .5e-5, 3))
        kt = 2
        zpkt = (zt, pt, kt)
        self.zr, self.pr, self.kr = ds.cancelPZ(zpkt, tol=1e-5)

    def test_zr(self):
        self.assertTrue(np.allclose(self.zr, (1.0, ), atol=1e-8, rtol=1e-6))

    def test_pr(self):
        self.assertTrue(
            np.allclose(self.pr, (1.0 + 2e-5, ), atol=1e-8, rtol=1e-6))

    def test_kr(self):
        self.assertTrue(np.allclose(self.kr, 2.0, atol=1e-8, rtol=1e-6))
