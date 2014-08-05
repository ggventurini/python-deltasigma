# -*- coding: utf-8 -*-
# test_bquantize.py
# Bipolar quantization test module
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

"""This is the test module for bquantize()"""

from __future__ import division, print_function
import unittest
import numpy as np
import deltasigma as ds
import scipy.io
import pkg_resources


class TestBQuantize(unittest.TestCase):
    """Test class for bquantize()"""
    def setUp(self):
        self.ran = np.arange(100)
    
    def test_bquantize(self):
        """Test function for bquantize()
        """
        x = np.linspace(-10, 10, 101)
        y = ds.bquantize(x)
        yval = [yi.val for yi in y]
        ycsd = [yi.csd for yi in y]
        fname = pkg_resources.resource_filename(__name__, "test_data/test_bquantize.mat")
        s = scipy.io.loadmat(fname)['s']
        mval = []
        mcsd = []
        for i in range(s.shape[1]):
            mval.append(float(s[0, i][0]))
            mcsd.append(s[0, i][1])
        for i in range(len(mval)):
            self.assertTrue(np.allclose(mval[i], yval[i], atol=1e-8, rtol=1e-5))
            self.assertTrue(np.prod(mcsd[i].shape) + np.prod(ycsd[i].shape) == 0 or \
                   mcsd[i].shape == ycsd[i].shape)
            if 0 not in ycsd[i].shape:
                self.assertTrue(np.allclose(mcsd[i], ycsd[i], atol=1e-8, rtol=1e-5))

