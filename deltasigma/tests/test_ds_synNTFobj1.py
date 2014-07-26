# -*- coding: utf-8 -*-
# test_ds_synNTFobj1.py
# This module provides the tests for the ds_synNTFobj1() function.
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


class TestDSSynNTFObj1(unittest.TestCase):
    """Test functions for ds_synNTFobj1()"""

    def setUp(self):
        pass

    def test_ds_synNTFobj1_1(self):
        """Test function for ds_synNTFobj1() 1/2"""
        res = -27.167735573627283
        tv = ds.ds_synNTFobj1(.5, (.9, 2), 64, .1)
        self.assertTrue(np.allclose((res,), (tv, ), atol=1e-8, rtol=1e-5))

    def test_ds_synNTFobj1_2(self):
        """Test function for ds_synNTFobj1() 2/2"""
        res = -43.0365
        tv = ds.ds_synNTFobj1(.5, (.9, 2), 64, 0.)
        self.assertTrue(np.allclose((res,), (tv, ), atol=1e-8, rtol=1e-5))

