# -*- coding: utf-8 -*-
# test_rmsGain.py
# This module provides the tests for the rmsGain function.
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

"""This module provides the test class for the rmsGain() function.
"""

import unittest
import numpy as np
import deltasigma as ds

class TestRmsGain(unittest.TestCase):
    """Test class for rmsGain()"""

    def setUp(self):
        num = (1,)
        den = (1, 2, 10)
        self.H = (num, den)
        self.f1 = 0.001
        self.f2 = 0.5
        self.res1 = 0.102245275091

    def test_rmsGain(self):
        """Test function for rmsGain()"""
        res = ds.rmsGain(self.H, self.f1, self.f2, N=1000)
        self.assertTrue(np.allclose((res,), (self.res1,), rtol=1e-05,
                                    atol=1e-08))

