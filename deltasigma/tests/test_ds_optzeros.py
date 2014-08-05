# -*- coding: utf-8 -*-
# test_ds_optzeros.py
# This module provides the tests for the ds_freq() function.
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
import scipy.io
from os.path import join


class TestDSOptZeros(unittest.TestCase):
    """Test function for ds_optzeros()"""

    def setUp(self):
        file_path = join('test_data', 'test_ds_optzeros.mat')
        fname = pkg_resources.resource_filename(__name__, file_path)
        self.res = scipy.io.loadmat(fname)['res']

    def test_opt_zeros(self):
        """ Test 14 different optzeros() calls. """
        ns = ('n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10',
              'n11', 'n12', 'n13', 'n14')

        for i, n in enumerate(ns):
            for opt in (0, 1, 2):
                self.assertTrue(
                    np.allclose(
                        self.res[n][0][0][:, opt], ds.ds_optzeros(i+1, opt),
                        atol=1e-10, rtol=1e-6))
