# -*- coding: utf-8 -*-
# test_calculateTF.py
# This module provides the tests for the calculateTF() function.
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
from deltasigma._utils import cplxpair


class TestCalculateTF(unittest.TestCase):
    """Test function for calculateTF()"""
    def setUp(self):
        ABCD = [[1.0, 0.0, 0.0, 0.044408783846879, -0.044408783846879],
                [0.999036450096481, 0.997109907515262, -0.005777399147297,
                 0.0, 0.499759089304780],
                [0.499759089304780, 0.999036450096481, 0.997109907515262,
                 0.0, -0.260002096136488],
                [0.0, 0.0, 1.0,  0.0, -0.796730400347216]]
        ABCD = np.array(ABCD)
        (ntf, stf) = ds.calculateTF(ABCD)
        (ntf_zeros, ntf_poles) = (np.roots(ntf.num), np.roots(ntf.den))
        (stf_zeros, stf_poles) = (np.roots(stf.num), np.roots(stf.den))
        mntf_poles = np.array((1.498975311463384, 1.102565142679772,
                               0.132677264750882))
        mntf_zeros = np.array((0.997109907515262 + 0.075972576202904j,
                               0.997109907515262 - 0.075972576202904j,
                               1.000000000000000 + 0.000000000000000j))
        mstf_zeros = np.array((-0.999999999999996,))
        mstf_poles = np.array((1.498975311463384, 1.102565142679772,
                               0.132677264750882))

        # for some reason, sometimes the zeros are in different order.
        (self.ntf_zeros, self.mntf_zeros) = (cplxpair(ntf_zeros),
                                             cplxpair(mntf_zeros))
        (self.stf_zeros, self.mstf_zeros) = (cplxpair(stf_zeros),
                                             cplxpair(mstf_zeros))
        (self.ntf_poles, self.mntf_poles) = (cplxpair(ntf_poles),
                                             cplxpair(mntf_poles))
        (self.stf_poles, self.mstf_poles) = (cplxpair(stf_poles),
                                             cplxpair(mstf_poles))

    def test_ntf_zeros_equal(self):
        self.assertTrue(
            np.allclose(self.ntf_zeros, self.mntf_zeros, rtol=1e-5, atol=1e-8))

    def test_ntf_poles_equal(self):
        self.assertTrue(
            np.allclose(self.ntf_poles, self.mntf_poles, rtol=1e-5, atol=1e-8))

    def test_stf_zeros_equal(self):
        self.assertTrue(
            np.allclose(self.stf_zeros, self.mstf_zeros, rtol=1e-5, atol=1e-8))

    def test_stf_poles_equal(self):
        self.assertTrue(
            np.allclose(self.stf_poles, self.mstf_poles, rtol=1e-5, atol=1e-8))
