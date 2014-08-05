# -*- coding: utf-8 -*-
# test_scaleABCD.py
# This module provides the tests for the scaleABCD function.
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

"""This module provides the test class for the scaleABCD() function.
"""

import unittest
import numpy as np
import deltasigma as ds

class TestScaleABCD(unittest.TestCase):
    """Test class for scaleABCD()"""

    def setUp(self):
        order = 8
        osr = 32
        self.nlev = 2
        self.f0 = 0.125
        Hinf = 1.5
        form = 'CRFB'
        # Optimized zero placement
        ntf = ds.synthesizeNTF(order, osr, 2, Hinf, self.f0)
        a, g, b, c = ds.realizeNTF(ntf, form)
        # we pass b = b[0]
        # if you use a single feed-in for the input, b may be scalar too
        ABCD = ds.stuffABCD(a, g, b[0], c, form)
        # now we correct b for the assert
        b = np.concatenate(( # Use a single feed-in for the input
                            np.atleast_1d(b[0]),
                            np.zeros((max(b.shape) - 1,))
                          ))
        self.ABCD0 = ABCD.copy()

        # References
        self.Sdiag_ref = np.array([71.9580, 51.9359, 8.2133, 6.5398,
                                   1.9446, 1.2070, 0.4223, 0.3040])
        self.umax_ref = 0.8667
        ABCD_ref1 = np.array([[1., -0.7320, 0, 0, 0, 0],
                              [0.7218, 0.4717, 0, 0, 0, 0],
                              [0, 0.1581, 1., -0.7357, 0, 0],
                              [0, 0.1259, 0.7962, 0.4142, 0, 0],
                              [0, 0, 0, 0.2973, 1., -0.9437],
                              [0, 0, 0, 0.1846, 0.6207, 0.4142],
                              [0, 0, 0, 0, 0, 0.3499],
                              [0, 0, 0, 0, 0, 0.2518],
                              [0, 0, 0, 0, 0, 0]])
        ABCD_ref2 = np.array([[0, 0, 0.0858, -0.0858],
                              [0, 0, 0.0619, 0.0428],
                              [0, 0, 0, 0.0642],
                              [0, 0, 0, 0.1835],
                              [0, 0, 0, 0.2447],
                              [0, 0, 0, 0.0581],
                              [1., -0.8971, 0, -0.0076],
                              [0.7197, 0.3543, 0, -0.1746],
                              [0, 3.29, 0, 0]])
        self.ABCD_ref = np.hstack((ABCD_ref1, ABCD_ref2))

    def test_scaleABCD(self):
        """Test function for scaleABCD()"""
        ABCD, umax, S = ds.scaleABCD(self.ABCD0, self.nlev, self.f0)
        # mapping the NTF to states, there is not a perfect match between 
        # the original code and the scipy version. -> rtol approx 20%
        #if not np.allclose(ABCD, ABCD_ref, atol=1e-2, rtol=3e-1):
        #    aerr = ABCD_ref-ABCD
        #    rerr = 2*(ABCD_ref-ABCD)/(ABCD_ref+ABCD)
        #    print(repr(ABCD_ref))
        #    print(repr(ABCD))
        #    print(aerr)
        #    print(rerr)
        # this is a rather high relative error. We get it on Travis-CI
        # Probably related to the libs used?
        self.assertTrue(np.allclose(ABCD, self.ABCD_ref, atol=1e-2, rtol=30e-2))
        self.assertTrue(np.allclose(umax, self.umax_ref, atol=1e-4, rtol=1e-3))
        self.assertTrue(np.allclose(np.diag(S), self.Sdiag_ref, atol=1e-2, rtol=25e-1))

        

