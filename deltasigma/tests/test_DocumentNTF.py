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
        # Note - nlev is never used
        # Was there another originally here where it would be used?
        nlev = 2
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

    def test_first_ntf_plot(self):
        """ Check plot with f0 = 0 """
        # check that DocumentNTF plots with no errors.
        f0 = 0.0
        self.assertIsNone(ds.DocumentNTF(self.ABCD, self.osr, f0))

    def test_second_ntf_plot(self):
        """ Check plot with f0 = 0.333 """
        f0 = 0.333
        self.assertIsNone(ds.DocumentNTF(self.ntf2, self.osr, f0))
