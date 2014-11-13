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

from deltasigma import synthesizeNTF, realizeNTF, stuffABCD, cplxpair

class TestSimulateDSM(unittest.TestCase):
    """Test class for simulateDSM()"""

    def setUp(self):
        fname = pkg_resources.resource_filename(__name__,
                "test_data/test_simulateDSM.mat")
        self.v_ref = scipy.io.loadmat(fname)['v']
        self.xn_ref = scipy.io.loadmat(fname)['xn']
        self.xmax_ref = scipy.io.loadmat(fname)['xmax']
        self.y_ref = scipy.io.loadmat(fname)['y']
        self.u_ref = scipy.io.loadmat(fname)['u']
        self.ABCD_ref = scipy.io.loadmat(fname)['ABCD']
        self.zeros = cplxpair(scipy.io.loadmat(fname)['zeros'])
        self.poles = cplxpair(scipy.io.loadmat(fname)['poles'])
        self.v_ref = self.v_ref.reshape(-1)
        self.y_ref = self.y_ref.reshape(-1)
        OSR = 32
        self.H = synthesizeNTF(5, OSR, 10)
        N = 8192
        f = 85
        self.u = 0.5*np.sin(2*np.pi*f/N*np.arange(N))
        a, g, b, c = realizeNTF(self.H, 'CRFB')
        self.ABCD = stuffABCD(a, g, b, c, form='CRFB')

    def test_precheck(self):
        """Test function for simulateDSM_* preliminary 1/1"""
        from deltasigma._utils import _get_zpk
        zeros, poles, _ = _get_zpk(self.H)
        assert np.allclose(cplxpair(zeros), self.zeros, atol=1e-6, rtol=1e-4)
        assert np.allclose(cplxpair(poles), self.poles, atol=1e-6, rtol=1e-4)
        assert np.allclose(self.ABCD, self.ABCD_ref, atol=1e-5, rtol=1e-3)
        assert np.allclose(self.u, self.u_ref, atol=1e-5, rtol=1e-3)

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
        v, _, _, y = ds._simulateDSM._simulateDSM_cblas(self.u, self.H, nlev=10,
                                                        x0=0.0, store_xn=True,
                                                        store_xmax=True,
                                                        store_y=True)
        v = v.reshape(-1)
        y = y.reshape(-1)
        self.assertTrue(np.allclose(v[:1200], self.v_ref[:1200], atol=1e-6,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(y[:1200], self.y_ref[:1200], atol=1e-1,
                                    rtol=1e-2))
        v, xn, xmax, y = ds._simulateDSM._simulateDSM_cblas(self.u, self.ABCD,
                                                 nlev=10, x0=0, store_xn=True,
                                                store_xmax=True, store_y=True)
        v = v.reshape(-1)
        y = y.reshape(-1)
        self.assertTrue(np.allclose(v[:1200], self.v_ref[:1200], atol=1e-6,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(y[:1200], self.y_ref[:1200], atol=1e-1,
                                    rtol=1e-3))
        self.assertTrue(np.allclose(xmax, self.xmax_ref, atol=1e-1, rtol=1e-1))
        self.assertTrue(np.allclose(xn[:, :1200], self.xn_ref[:, :1200],
                                    atol=1e-1, rtol=1e-3))

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
        v, xn, xmax, y = ds._simulateDSM._simulateDSM_scipy_blas(self.u, self.H,
                                                  nlev=10, x0=0., store_xn=True,
                                                  store_xmax=True, store_y=True)
        v = v.reshape(-1)
        y = y.reshape(-1)
        self.assertTrue(np.allclose(v[:1200], self.v_ref[:1200], atol=1e-6,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(y[:1200], self.y_ref[:1200], atol=1e-1,
                                    rtol=1e-2))
        v, xn, xmax, y = ds._simulateDSM._simulateDSM_scipy_blas(self.u,
                                    self.ABCD, nlev=10, x0=0, store_xn=True,
                                    store_xmax=True, store_y=True)
        v = v.reshape(-1)
        y = y.reshape(-1)
        self.assertTrue(np.allclose(v[:1200], self.v_ref[:1200], atol=1e-6,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(y[:1200], self.y_ref[:1200], atol=1e-1,
                                    rtol=1e-3))
        self.assertTrue(np.allclose(xmax, self.xmax_ref, atol=1e-1, rtol=1e-1))
        self.assertTrue(np.allclose(xn[:, :1200], self.xn_ref[:, :1200],
                                    atol=1e-1, rtol=1e-3))

    def test_simulateDSM_python(self):
        """Test function for simulateDSM_python()"""
        v, xn, xmax, y = ds._simulateDSM._simulateDSM_python(self.u, self.H,
                                                             nlev=10., x0=0.)
        v = v.reshape(-1)
        y = y.reshape(-1)
        self.assertTrue(np.allclose(v[:1200], self.v_ref[:1200], atol=1e-6,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(y[:1200], self.y_ref[:1200], atol=1e-1,
                                    rtol=1e-2))
        v, xn, xmax, y = ds._simulateDSM._simulateDSM_python(self.u, self.ABCD,
                                                             nlev=10, x0=0.)
        v = v.reshape(-1)
        y = y.reshape(-1)
        self.assertTrue(np.allclose(v[:1200], self.v_ref[:1200], atol=1e-6,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(y[:1200], self.y_ref[:1200], atol=1e-1,
                                    rtol=1e-3))
        self.assertTrue(np.allclose(xmax, self.xmax_ref, atol=1e-1, rtol=1e-1))
        self.assertTrue(np.allclose(xn[:, :1200], self.xn_ref[:, :1200],
                                    atol=1e-1, rtol=1e-3))

    def test_simulateDSM_scipy_cblas(self):
        """Test function for scipy vs cblas"""
        # skip if not available
        try:
            from nose.plugins.skip import SkipTest
        except ImportError:
            SkipTest = None
        if ds._simulateDSM._simulateDSM_scipy_blas is None or \
            ds._simulateDSM._simulateDSM_cblas is None:
            if SkipTest is not None:
                raise SkipTest
            return
        vb, _, _, yb = ds._simulateDSM._simulateDSM_scipy_blas(self.u,
                                    self.H, nlev=10, x0=0.,
                                    store_xn=True, store_xmax=True,
                                    store_y=True)
        vc, _, _, yc = ds._simulateDSM._simulateDSM_cblas(self.u, self.H,
                                    nlev=10, x0=0., store_xn=True,
                                    store_xmax=True, store_y=True)
        self.assertTrue(np.allclose(vb, vc, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(yb, yc, atol=1e-6, rtol=1e-4))
        vb, xnb, xmaxb, yb = ds._simulateDSM._simulateDSM_scipy_blas(self.u,
                                    self.ABCD, nlev=10, x0=0., store_xn=True,
                                    store_xmax=True, store_y=True)
        vc, xnc, xmaxc, yc = ds._simulateDSM._simulateDSM_cblas(self.u,
                                    self.ABCD, nlev=10, x0=0., store_xn=True,
                                    store_xmax=True, store_y=True)
        self.assertTrue(np.allclose(vb, vc, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(xnb, xnc, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(xmaxb, xmaxc, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(yb, yc, atol=1e-6, rtol=1e-4))

    def test_simulateDSM_python_scipy(self):
        """Test function for Python vs scipy"""
        # skip if not available
        try:
            from nose.plugins.skip import SkipTest
        except ImportError:
            SkipTest = None
        if ds._simulateDSM._simulateDSM_scipy_blas is None:
            if SkipTest is not None:
                raise SkipTest
            return
        vb, _, _, yb = ds._simulateDSM._simulateDSM_scipy_blas(self.u,
                                    self.H, nlev=10, x0=0.,
                                    store_xn=True, store_xmax=True,
                                    store_y=True)
        vc, _, _, yc = ds._simulateDSM._simulateDSM_python(self.u, self.H,
                                    nlev=10, x0=0.)
        self.assertTrue(np.allclose(vb, vc, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(yb, yc, atol=1e-6, rtol=1e-4))
        vb, xnb, xmaxb, yb = ds._simulateDSM._simulateDSM_scipy_blas(self.u,
                                    self.ABCD, nlev=10, x0=0., store_xn=True,
                                    store_xmax=True, store_y=True)
        vc, xnc, xmaxc, yc = ds._simulateDSM._simulateDSM_python(self.u,
                                    self.ABCD, nlev=10, x0=0.)
        self.assertTrue(np.allclose(vb, vc, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(xnb, xnc, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(xmaxb, xmaxc, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(yb, yc, atol=1e-6, rtol=1e-4))

    def test_simulateDSM_python_cblas(self):
        """Test function for Python vs cblas"""
        # skip if not available
        try:
            from nose.plugins.skip import SkipTest
        except ImportError:
            SkipTest = None
        if ds._simulateDSM._simulateDSM_cblas is None:
            if SkipTest is not None:
                raise SkipTest
            return
        vb, _, _, yb = ds._simulateDSM._simulateDSM_cblas(self.u,
                                    self.H, nlev=10, x0=0.,
                                    store_xn=True, store_xmax=True,
                                    store_y=True)
        vc, _, _, yc = ds._simulateDSM._simulateDSM_python(self.u, self.H,
                                    nlev=10, x0=0.)
        self.assertTrue(np.allclose(vb, vc, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(yb, yc, atol=1e-6, rtol=1e-4))
        vb, xnb, xmaxb, yb = ds._simulateDSM._simulateDSM_cblas(self.u,
                                    self.ABCD, nlev=10, x0=0., store_xn=True,
                                    store_xmax=True, store_y=True)
        vc, xnc, xmaxc, yc = ds._simulateDSM._simulateDSM_python(self.u,
                                    self.ABCD, nlev=10, x0=0.)
        self.assertTrue(np.allclose(vb, vc, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(xnb, xnc, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(xmaxb, xmaxc, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(yb, yc, atol=1e-6, rtol=1e-4))

