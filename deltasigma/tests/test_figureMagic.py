# -*- coding: utf-8 -*-
# test_figureMagic.py
# This module provides the tests for the figureMagic function.
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
import pylab as plt


class TestFigureMagic(unittest.TestCase):
    """Test functions for figureMagic()"""

    def setUp(self):
        pass

    def test_figureMagic(self):
        """test plotting - None should be returned."""
        a = np.arange(10)
        plt.figure()
        plt.plot(a)
        self.assertIsNone(
            ds.figureMagic(
                xRange=[1, 10], dx=1, xLab=None, yRange=[2, 8], dy=.5,
                yLab=None, size=(10, 6), name="Test plot"))
