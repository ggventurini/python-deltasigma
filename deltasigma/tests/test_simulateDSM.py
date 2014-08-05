# -*- coding: utf-8 -*-
# test_simulateDSM.py
# This module provides the tests for the simulateDSM function.
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

"""This module provides the test class for the simulateDSM() function.
"""

import unittest
import pkg_resources

import numpy as np
import deltasigma as ds
import scipy.io

from deltasigma import synthesizeNTF, realizeNTF, stuffABCD

class TestSimulateDSM(unittest.TestCase):
    """Test class for simulateDSM()"""

    def setUp(self):
        fname = pkg_resources.resource_filename(__name__, "test_data/test_simulateDSM.mat")
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

    def test_simulateDSM_cblas(self):
        """Test function for simulateDSM_cblas()"""
        # skip if not available
        try:
            from nose.plugins.skip import SkipTest
        except ImportError:
            SkipTest = None
        if ds._simulateDSM._simulateDSM_cblas is None:
            if SkipTest is not None:
                raise SkipTest
            return
        v, xn, xmax, y = ds._simulateDSM._simulateDSM_cblas(self.u, self.H, nlev=2,
                                    x0=0, store_xn=True, store_xmax=True, store_y=True)
        self.assertTrue(np.allclose(v.reshape(-1), self.v_ref.reshape(-1),
                                    atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(y, self.y_ref, atol=1e-6, rtol=1e-4))
        v, xn, xmax, y = ds._simulateDSM._simulateDSM_cblas(self.u, self.ABCD, nlev=2,
                                    x0=0, store_xn=True, store_xmax=True, store_y=True)
        self.assertTrue(np.allclose(v, self.v_ref, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(y, self.y_ref, atol=1e-6, rtol=1e-4))

    def test_simulateDSM_scipy_blas(self):
        """Test function for simulateDSM_scipy_blas()"""
        # skip if not available
        try:
            from nose.plugins.skip import SkipTest
        except ImportError:
            SkipTest = None
        if ds._simulateDSM._simulateDSM_scipy_blas is None:
            if SkipTest is not None:
                raise SkipTest
            return
        v, xn, xmax, y = ds._simulateDSM._simulateDSM_scipy_blas(self.u, self.H, nlev=2,
                                    x0=0, store_xn=True, store_xmax=True, store_y=True)
        self.assertTrue(np.allclose(v.reshape(-1), self.v_ref.reshape(-1), atol=1e-6,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(y, self.y_ref, atol=1e-6, rtol=1e-4))
        v, xn, xmax, y = ds._simulateDSM._simulateDSM_scipy_blas(self.u, self.ABCD, nlev=2,
                                    x0=0, store_xn=True, store_xmax=True, store_y=True)
        self.assertTrue(np.allclose(v, self.v_ref, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(y, self.y_ref, atol=1e-6, rtol=1e-4))

    def test_simulateDSM_python(self):
        """Test function for simulateDSM_python()"""
        v, xn, xmax, y = ds._simulateDSM._simulateDSM_python(self.u, self.H, nlev=2, x0=0)
        self.assertTrue(np.allclose(v.reshape(-1), self.v_ref.reshape(-1), atol=1e-6,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(y, self.y_ref, atol=1e-6, rtol=1e-4))
        v, xn, xmax, y = ds._simulateDSM._simulateDSM_python(self.u, self.ABCD, nlev=2, x0=0)
        self.assertTrue(np.allclose(v, self.v_ref, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(y, self.y_ref, atol=1e-6, rtol=1e-4))

