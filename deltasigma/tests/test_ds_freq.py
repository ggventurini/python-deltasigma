# -*- coding: utf-8 -*-
# test_ds_freq.py
# This module provides the tests for the ds_freq() function.
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


class TestDSFreq(unittest.TestCase):
    """Tests ds_freq()"""

    def setUp(self):
        pass

    def test_res_equals_tres(self):
        """Check ds_freq() output to known values."""
        a = ds.ds_freq(osr=128, f0=0.0, quadrature=True)
        b = np.diff(a)
        res = (0.00190595677588, 0.00510204081633,
               0.207803148686, 0.00491921819577)
        tres = (a.mean(), b.mean(), a.std(), b.std())
        self.assertTrue(np.allclose(res, tres, atol=1e-8, rtol=1e-5))
