# -*- coding: utf-8 -*-
# test_simulateDSM.py
# This module provides the tests for the simulateDSM_numba function.
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

"""This module provides the test class for the simulateDSM() function.
"""

import unittest
import pkg_resources

import numpy as np
import deltasigma as ds
import scipy.io

from deltasigma import synthesizeNTF, realizeNTF, stuffABCD

class TestSimulateDSM(unittest.TestCase):
    """Test class for simulateDSM_numba()"""

    def setUp(self):
        fname = pkg_resources.resource_filename(
            __name__, "test_data/test_simulateDSM.mat")
        self.v_ref = scipy.io.loadmat(fname)['v']
        self.xn_ref = scipy.io.loadmat(fname)['xn']
        self.xmax_ref = scipy.io.loadmat(fname)['xmax']
        self.y_ref = scipy.io.loadmat(fname)['y']
        OSR = 32
        self.H = synthesizeNTF(5, OSR, 1)
        N = 8192
        f = 85
        self.u = 0.5*np.sin(2*np.pi*f/N*np.arange(N))
        a, g, b, c = realizeNTF(self.H, 'CRFB')
        self.ABCD = stuffABCD(a, g, b, c, form='CRFB')

    def test_simulateDSM_numba_ABCD(self):
        """Test function for simulateDSM_numba()"""
        #print self.ABCD
        #print self.u
        #print self.H
        v, xn, xmax, y = ds._simulateDSM._simulateDSM_numba(
            self.u, self.ABCD, 2, 0)
        self.assertTrue(np.allclose(v, self.v_ref, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(y, self.y_ref, atol=1e-6, rtol=1e-4))

#    def test_simulateDSM_numba_H(self):
#        v, xn, xmax, y = ds._simulateDSM._simulateDSM_numba(
#            self.u, self.H, 2, 0)
#        self.assertTrue(np.allclose(
#            v.reshape(-1), self.v_ref.reshape(-1), atol=1e-6, rtol=1e-4))
#        self.assertTrue(np.allclose(y, self.y_ref, atol=1e-6, rtol=1e-4))
