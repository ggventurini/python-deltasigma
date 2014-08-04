# -*- coding: utf-8 -*-
# test_SIunits.py
# This module provides the tests for the SIunits function.
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

"""This module provides the test class for the SIunits() function.
"""

import unittest
import numpy as np

from deltasigma import SIunits

class TestSIunits(unittest.TestCase):
    """Test class for SIunits()"""

    def setUp(self):
        pass

    def test_SIunits(self):
        """Test function for SIunits() 1/2"""
        tv = (0, 1, 1e3, 2100312.24, .32545, 21e-9, 34e-12, 9569300e-12)
        correct = (0, ''), (1, ''), (1e3, 'k'), (1e6, 'M'), (1e-3, 'm'), \
                  (1e-9, 'n'), (1e-12, 'p'), (1e-6, 'u')
        f, p = SIunits(tv)
        res = zip(f,p)
        for r, c in zip(res, correct):
            self.assertTrue(r[0] == c[0] and r[1] == c[1])

    def test_SIunits_2(self):
        """Test function for SIunits() 2/2"""
        # test scalars
        tv = 2100312.24
        correct = (1e6, 'M')
        f, p = SIunits(tv)
        self.assertTrue(f == correct[0] and p == correct[1])

