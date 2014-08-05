# -*- coding: utf-8 -*-
# test_circ_smooth.py
# This module provides the tests for the circ_smooth() function.
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
import pkg_resources
from scipy.io import loadmat
from os.path import join


class TestCircSmooth(unittest.TestCase):
    """Test function for circ_smooth()"""

    def setUp(self):
        file_path = join('test_data', 'test_circ_smooth.mat')
        fname = pkg_resources.resource_filename(__name__, file_path)
        self.bt = loadmat(fname)['b']

    def test_circ_smooth(self):
        A = np.arange(1, 101)
        b = ds.circ_smooth(A, 16)
        self.assertTrue(np.allclose(self.bt, b, atol=1e-8, rtol=1e-5))
