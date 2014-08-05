
# -*- coding: utf-8 -*-
# test_undbp.py
# This module provides the tests for the undbp function.
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

"""This module provides the test class for the undbp() function.
"""

import unittest
import numpy as np

from deltasigma import undbp

class TestUndbp(unittest.TestCase):
    """Test class for undbp()"""

    def setUp(self):
        pass

    def test_undbp_1(self):
        """Test function for undbp() 1/2"""
        self.assertTrue(np.allclose([undbp(53.05)], [201836.636368], rtol=1e-05,
                                    atol=1e-08))

    def test_undbp_2(self):
        """Test function for undbp() 2/2"""
        self.assertTrue(np.allclose([undbp(3)], [1.99526231497], rtol=1e-05,
                                    atol=1e-08))

