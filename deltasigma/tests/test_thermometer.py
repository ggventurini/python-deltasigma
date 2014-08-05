# -*- coding: utf-8 -*-
# test_thermometer.py
# This module provides the tests for the thermometer function.
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

"""This module provides the test class for the thermometer() function.
"""

import unittest
import numpy as np

from deltasigma import thermometer

class TestThermometer(unittest.TestCase):
    """Test class for thermometer()"""

    def setUp(self):
        self.tv = np.arange(50)
        self.rm = np.zeros((70, 50))
        for i in range(50):
            self.rm[:i, i] = np.ones(self.rm[:i, i].shape)

    def test_thermometer(self):
        """Test function for thermometer() 1/3"""
        self.assertTrue(np.allclose(thermometer(self.tv, 70), self.rm,
                                    rtol=1e-05, atol=1e-08))

