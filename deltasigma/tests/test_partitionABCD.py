# -*- coding: utf-8 -*-
# test_partitionABCD.py
# This module provides the tests for the partitionABCD function.
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

"""This module provides the test class for the partitionABCD() function.
"""

import unittest
import numpy as np
import deltasigma as ds

from nose.tools import raises
from scipy.signal import lti

class TestPartitionABCD(unittest.TestCase):
    """Test class for partitionABCD()"""

    def setUp(self):
        # data for test 1
        self.ob = lti((1, ), (1, 2, 10)).to_ss()
        ab = np.hstack((self.ob.A, self.ob.B))
        cd = np.hstack((self.ob.C, self.ob.D.reshape((1,1))))
        self.ABCD1 = np.vstack((ab, cd))
        # data for test2
        ABCD = [[1.000000000000000, 0., 0., 0.044408783846879, -0.044408783846879],
                [0.999036450096481, 0.997109907515262, -0.005777399147297, 0., 0.499759089304780],
                [0.499759089304780, 0.999036450096481, 0.997109907515262,  0., -0.260002096136488],
                [0, 0, 1., 0, -0.796730400347216]]
        self.ABCD2 = np.array(ABCD)
        self.at2 = self.ABCD2[:3, :3]
        self.bt2 = self.ABCD2[:3, 3:]
        self.ct2 = self.ABCD2[3, :3]
        self.dt2 = self.ABCD2[3, 3:]

        self.at3 = self.ABCD2[:2, :2]
        self.bt3 = self.ABCD2[:2, 2:]
        self.ct3 = self.ABCD2[2:, :2]
        self.dt3 = self.ABCD2[2:, 2:]

    def test_partitionABCD_1(self):
        """Test function for partitionABCD() 1/5"""
        a, b, c, d = ds.partitionABCD(self.ABCD1)
        self.assertTrue(np.allclose(a, self.ob.A, rtol=1e-5, atol=1e-8))
        self.assertTrue(np.allclose(b, self.ob.B, rtol=1e-5, atol=1e-8))
        self.assertTrue(np.allclose(c, self.ob.C, rtol=1e-5, atol=1e-8))
        self.assertTrue(np.allclose(d, self.ob.D, rtol=1e-5, atol=1e-8))

    def test_partitionABCD_2(self):
        """Test function for partitionABCD() 2/5"""
        ar, br, cr, dr = ds.partitionABCD(self.ABCD2, m=2)
        self.assertTrue(np.allclose(self.at2, ar, rtol=1e-5, atol=1e-8))
        self.assertTrue(np.allclose(self.bt2, br, rtol=1e-5, atol=1e-8))
        self.assertTrue(np.allclose(self.ct2, cr, rtol=1e-5, atol=1e-8))
        self.assertTrue(np.allclose(self.dt2, dr, rtol=1e-5, atol=1e-8))

    def test_partitionABCD_3(self):
        """Test function for partitionABCD() 3/5"""
        ar, br, cr, dr = ds.partitionABCD(self.ABCD2, m=2, r=1)
        self.assertTrue(np.allclose(self.at2, ar, rtol=1e-5, atol=1e-8))
        self.assertTrue(np.allclose(self.bt2, br, rtol=1e-5, atol=1e-8))
        self.assertTrue(np.allclose(self.ct2, cr, rtol=1e-5, atol=1e-8))
        self.assertTrue(np.allclose(self.dt2, dr, rtol=1e-5, atol=1e-8))

    def test_partitionABCD_4(self):
        """Test function for partitionABCD() 4/5"""
        ar, br, cr, dr = ds.partitionABCD(self.ABCD2, m=3, r=2)
        self.assertTrue(np.allclose(self.at3, ar, rtol=1e-5, atol=1e-8))
        self.assertTrue(np.allclose(self.bt3, br, rtol=1e-5, atol=1e-8))
        self.assertTrue(np.allclose(self.ct3, cr, rtol=1e-5, atol=1e-8))
        self.assertTrue(np.allclose(self.dt3, dr, rtol=1e-5, atol=1e-8))

    @raises(ValueError)
    def test_partitionABCD_5(self):
        """Test function for partitionABCD() 5/5"""
        ar, br, cr, dr = ds.partitionABCD(self.ABCD2, m=2, r=2)

