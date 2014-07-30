# -*- coding: utf-8 -*-
# test_evalRPoly.py
# This module provides the tests for the evalRPoly function.
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

class TestEvalRPoly(unittest.TestCase):
    """Test function for evalRPoly()"""
    def setUp(self):
        pass

    def test_evalRPoly(self):
        """test_evalRPoly"""
        x = np.arange(1001) - 500
        a = [1, 0, 1, 2]
        r1 = np.polyval(a, x)
        rts = np.roots(a)
        r2 = ds.evalRPoly(rts, x)
        self.assertTrue(np.allclose(r1, r2, atol=1e-8, rtol=1e-5))
