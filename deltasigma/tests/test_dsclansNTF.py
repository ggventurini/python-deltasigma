# -*- coding: utf-8 -*-
# test_dsclansNTF.py
# This module provides the tests for dsclansNTF() function.
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


class TestDSClansNTF(unittest.TestCase):
    """Test functions for dsclansNTF()"""

    def setUp(self):
        self.x1 = ds.dsclansNTF(np.arange(1, 100.001, .001), 3, .5, 100)

    def test_rp_rt(self):
        """ Test rp == rt  """
        rt = np.array(self.x1[1][:])
        rt.sort()
        rp = np.array([0, -0.016805373426715, 0.014809370415763])
        rp.sort()
        self.assertTrue(np.allclose(rp, rt, rtol=1e-5, atol=1e-6))

    def test_rz_x0(self):
        """ Test rz == x[0] """
        rz = (100., )
        self.assertTrue(np.allclose(rz, self.x1[0], rtol=1e-5, atol=1e-8))

    def test_rp_tp(self):
        """ Test rp == tp """
        rp = np.array(
            [0.35378443 + 0.0j, 0.34187718 + 0.22831019j,
             0.34187718 - 0.22831019j, 0.32978826 + 0.59355161j,
             0.32978826 - 0.59355161j])
        x = np.array([0.67623674, 0.9277613, 0.70365961, 0.60374009,
                      0.78008118])
        Hz = np.array([0.99604531 - 0.08884669j, 0.99604531 + 0.08884669j,
                       0.99860302 - 0.05283948j, 0.99860302 + 0.05283948j,
                       1.00000000 + 0.0j])
        tp = ds.dsclansNTF(x, order=5, rmax=.95, Hz=Hz)[1]
        self.assertTrue(np.allclose(rp, tp, rtol=1e-5, atol=1e-6))
