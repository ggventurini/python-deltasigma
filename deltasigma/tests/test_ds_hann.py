# -*- coding: utf-8 -*-
# test_ds_hann.py
# This module provides the tests for the ds_hann() function.
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

"""This module provides the test class for the ds_hann() function.
"""

import unittest
import numpy as np
import deltasigma as ds

class TestDs_Hann(unittest.TestCase):
    """Test class for ds_hann()"""
    def setUp(self):
        self.res = np.array([0.        ,  0.02148628,  0.06768441,  0.09549150,
                             0.06533781, -0.03015369, -0.1545085 , -0.24133259,
                             -0.22851372, -0.0954915 ])

    def test_ds_hann(self):
        """Test function for ds_hann()"""
        self.assertTrue(np.allclose(self.res, np.hanning(10) - ds.ds_hann(10), atol=1e-8, rtol=1e-5))
        

