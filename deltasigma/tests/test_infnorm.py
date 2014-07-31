# -*- coding: utf-8 -*-
# test_infnorm.py
# This module provides the tests for the infnorm function.
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


class TestInfNorm(unittest.TestCase):
    """Test functions for figureMagic()"""
    def setUp(self):
        num, den = np.poly([3, 0.3, 1]), np.poly([2, 0.5, .25])
        H = (num, den)
        self.Hinf, self.fmax = ds.infnorm(H)

    def test_infnorm_Hinf(self):
        """Test function for infnorm() checking Hinf"""
        self.assertTrue(np.allclose(
            self.Hinf, 1.84888889, atol=1e-8, rtol=1e-5))

    def test_infnorm_fmax(self):
        """Test function for infnorm() checking fmax"""
        self.assertTrue(np.allclose(
            self.fmax, 3.141592653589793/2.0/np.pi, atol=1e-8, rtol=1e-5))
