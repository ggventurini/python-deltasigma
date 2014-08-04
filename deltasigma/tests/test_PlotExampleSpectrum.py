# -*- coding: utf-8 -*-
# test_PlotExampleSpectrum.py
# This module provides the tests for the PlotExampleSpectrum function.
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

"""This module provides the test class for the PlotExampleSpectrum() function.
"""

import unittest
import numpy as np
import deltasigma as ds

class TestPlotExampleSpectrum(unittest.TestCase):
    """Test class for PlotExampleSpectrum()"""

    def setUp(self):
        pass

    def test_PlotExampleSpectrum(self):
        """Test function for PlotExampleSpectrum()"""
        order = 3
        osr = 32
        f0 = 0.
        Hinf = 1.5
        ntf = ds.synthesizeNTF(order, osr, 0, Hinf, f0)
        ret = ds.PlotExampleSpectrum(ntf, M=1, osr=osr, f0=f0)
        self.assertIsNone(ret)

