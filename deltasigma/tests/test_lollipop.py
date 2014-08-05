# -*- coding: utf-8 -*-
# test_lollipop.py
# This module provides the tests for the lollipop function.
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

"""This module provides the test class for the lollipop() function.
"""

import unittest
import numpy as np
import pylab as plt

from warnings import catch_warnings

from deltasigma import lollipop

class TestLollipop(unittest.TestCase):
    """Test class for lollipop()"""

    def setUp(self):
        pass

    def test_lollipop(self):
        """Test function for lollipop()"""
        t = np.arange(1, 20)*1e-3
        f = 20.
        a = np.sin(2*np.pi*f*t)
        plt.figure()
        with catch_warnings(record=True) as w:
            lollipop(t, a, color=None, lw=1.5, ybot=0.1)
            self.assertTrue(len(w) > 0)

