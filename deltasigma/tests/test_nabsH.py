# -*- coding: utf-8 -*-
# test_nabsH.py
# This module provides the tests for the nabsH function.
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

"""This module provides the test class for the nabsH() function.
"""

import unittest
import numpy as np
import deltasigma as ds

from nose.tools import raises

class TestNabsH(unittest.TestCase):
    """Test class for nabsH()"""

    def setUp(self):
        pass

    def test_nabsH(self):
        """Test function for nabsH()"""
        H = ([1, 2], [2, 0, .25], 1)
        N = 129
        w = np.linspace(0, 2*np.pi, num=N, endpoint=True)
        z = np.exp(1j*w)
        r1 = -np.abs(ds.evalTF(H, z))
        r2 = ds.nabsH(w, H)
        self.assertTrue(np.allclose(r1, r2, atol=1e-8, rtol=1e-5))

