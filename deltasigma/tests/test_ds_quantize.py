# -*- coding: utf-8 -*-
# test_ds_quantize.py
# This module provides the tests for the ds_quantize() function.
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


class TestDSQuantize(unittest.TestCase):
    """Test class for ds_quantize()"""

    def setUp(self):
        t = np.arange(-3, 3, .2)
        t = t.reshape((1, t.shape[0]))
        y = t
        for _ in range(2):
            y = np.concatenate((y, t), axis=0)
        n = np.array([2, 3, 4])
        self.re1_first = ds.ds_quantize(y, n)
        self.re1_second = ds.ds_quantize(y[1, :], n=3)
        self.re2 = np.array(
                [[-1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1.,
                  -1., -1., -1., -1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                   1., 1., 1., 1., 1., 1.],
                 [-2., -2., -2., -2., -2., -2., -2., -2., -2., -2., 0., 0.,
                   0., 0., 0., 0., 0., 0., 0., 0., 2., 2., 2., 2., 2., 2.,
                   2., 2., 2., 2.],
                 [-3., -3., -3., -3., -3., -1., -1., -1., -1., -1., -1., -1.,
                  -1., -1., -1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 3.,
                   3., 3., 3., 3.]])

    def test_first_quantize(self):
        """Test function for ds_quantize() 1/2 """
        self.assertTrue(
            np.allclose(self.re1_first, self.re2, atol=1e-8, rtol=1e-5))

    def test_second_quantizer(self):
        """Test function for ds_quantize() 2/2 """
        self.assertTrue(
            np.allclose(self.re1_second, self.re2[1, :], atol=1e-8, rtol=1e-5))
