# -*- coding: utf-8 -*-
# test_clans.py
# This module provides the tests for the clans() function.
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


class TestClans(unittest.TestCase):
    """Class doc string"""
    def setUp(self):
        self.ntf = ds.clans(5, 32, 5, .95, 1)
        self.poles = np.array((0.41835234+0.0j, 0.48922229+0.1709716j,
                               0.48922229-0.1709716j, 0.65244885+0.3817224j,
                               0.65244885-0.3817224j))

    def test_clans(self):
        """Test function for clans()"""
        self.assertTrue(
            np.allclose(self.poles, self.ntf[1], atol=1e-8, rtol=1e-5))
