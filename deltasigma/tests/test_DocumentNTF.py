# -*- coding: utf-8 -*-
# test_DocumentNTF.py
# This module provides the tests for the DocumentNTF() function.
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
from deltasigma._synthesizeNTF import synthesizeNTF
from deltasigma._realizeNTF import realizeNTF
from deltasigma._stuffABCD import stuffABCD


class TestDocumentNTF(unittest.TestCase):
    """Test function for DocumentNTF()"""

    def setUp(self):
        order = 4
        self.osr = 64
        f0 = 0.0
        Hinf = 1.5
        form = 'CRFB'
        ntf1 = synthesizeNTF(order, self.osr, 2, Hinf, f0)
        self.ntf2 = synthesizeNTF(order, self.osr, 2, Hinf, f0=.333)
        a, g, b, c = realizeNTF(ntf1, form)
        # Use a single feed-in for the input
        b = np.concatenate((
            np.atleast_1d(b[0]), np.zeros((max(b.shape) - 1))
            ))
        self.ABCD = stuffABCD(a, g, b, c, form)

    def test_documentNTF1(self):
        """Test function for DocumentNTF(): check plot with f0 = 0 1/2"""
        # check that DocumentNTF plots with no errors.
        f0 = 0.0
        self.assertIsNone(ds.DocumentNTF(self.ABCD, self.osr, f0))

    def test_documentNTF2(self):
        """Test function for DocumentNTF(): check plot with f0 = 0.333 2/2"""
        f0 = 0.333
        self.assertIsNone(ds.DocumentNTF(self.ntf2, self.osr, f0))

