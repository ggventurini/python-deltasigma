# -*- coding: utf-8 -*-
# test_plotPZ.py
# This module provides the tests for the plotPZ function.
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

"""This module provides the test class for the plotPZ() function.
"""

import unittest
import numpy as np
import pylab as plt

from deltasigma import plotPZ

class TestPlotPZ(unittest.TestCase):
    """Test class for plotPZ()"""

    def setUp(self):
        pass

    def test_plotPZ_1(self):
        """Test function for plotPZ() 1/2"""
        plt.figure()
        self.assertIsNone(plotPZ(((1, .2), (1, 0, .10)), color=('r', 'b'),
                                 showlist=True))

    def test_plotPZ_2(self):
        """Test function for plotPZ() 2/2"""
        plt.figure()
        self.assertIsNone(plotPZ(((1, 0, .10), (1, .2, .01)), showlist=False))

