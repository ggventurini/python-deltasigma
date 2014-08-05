# -*- coding: utf-8 -*-
# test_mapRtoQ.py
# This module provides the tests for the mapRtoQ function.
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

"""This module provides the test class for the mapRtoQ() function.
"""

from __future__ import division

import unittest
import numpy as np
import deltasigma as ds

class TestMapRtoQ(unittest.TestCase):
    """Test class for mapRtoQ()"""

    def setUp(self):
        self.test_matrix = np.arange(1, 25).reshape((4, 6)).T
        self.dq = np.array([[4.5 - 2.5j, 16.5 - 2.5j],
                            [6.5 - 2.5j, 18.5 - 2.5j],
                            [8.5 - 2.5j, 20.5 - 2.5j]])
        self.dp = np.array([[-3.5 + 4.5j, -3.5 + 16.5j],
                            [-3.5 + 6.5j, -3.5 + 18.5j],
                            [-3.5 + 8.5j, -3.5 + 20.5j]])

    def test_mapRtoQ(self):
        """Test function for mapRtoQ()"""
        resq, resp = ds.mapRtoQ(self.test_matrix)
        self.assertTrue(np.allclose(resq, self.dq, atol=1e-8, rtol=1e-5))
        self.assertTrue(np.allclose(resp, self.dp, atol=1e-8, rtol=1e-5))

