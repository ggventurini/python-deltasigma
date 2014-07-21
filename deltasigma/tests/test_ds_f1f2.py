# -*- coding: utf-8 -*-
# test_ds_f1f2.py
# This module provides the tests for the ds_f1f2 function.
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

"""This module provides the test class for the ds_f1f2() function.
"""

import unittest
import deltasigma as ds

class TestDs_f1f2(unittest.TestCase):
    """Test class for ds_f1f2()"""
    def setUp(self):
        self.f0 = 1e3
        self.OSR = 128
        self.cf = False

    def test_ds_f1f2_1(self):
        """Test function for ds_f1f2() 1/3"""
        t1f1, t1f2 = ds.ds_f1f2(self.OSR, self.f0, self.cf)
        self.assertTrue((t1f1, t1f2) == (self.f0 - 0.25/self.OSR, self.f0 + 0.25/self.OSR))

    def test_ds_f1f2_2(self):
        """Test function for ds_f1f2() 2/3"""
        t2f1, t2f2 = ds.ds_f1f2(self.OSR, self.f0, not self.cf)
        self.assertTrue((t2f1, t2f2) == (self.f0 - 0.5/self.OSR, self.f0 + 0.5/self.OSR))

    def test_ds_f1f2_3(self):
        """Test function for ds_f1f2() 3/3"""
        t3f1, t3f2 = ds.ds_f1f2(self.OSR, .24/self.OSR, self.cf)
        self.assertTrue((t3f1, t3f2) == (0., 0.5/self.OSR))
