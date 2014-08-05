# -*- coding: utf-8 -*-
# test_simulateSNR.py
# This module provides the tests for the simulateSNR function.
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

"""This module provides the test class for the simulateSNR() function.
"""

import unittest
import pkg_resources
import scipy.io
import numpy as np
import deltasigma as ds

from scipy.signal import lti

class TestSimulateSNR(unittest.TestCase):
    """Test class for simulateSNR()"""

    def setUp(self):
        pass

    def test_simulateSNR_1(self):
        """Test function for simulateSNR() 1/3"""
        # first test: f0 = 0
        # Load test references
        fname = pkg_resources.resource_filename(__name__,
                                                "test_data/test_snr_amp.mat")
        amp_ref = scipy.io.loadmat(fname)['amp'].reshape((-1,))
        snr_ref = scipy.io.loadmat(fname)['snr'].reshape((-1,))
        amp_user_ref = scipy.io.loadmat(fname)['amp_user'].reshape((-1,))
        snr_user_ref = scipy.io.loadmat(fname)['snr_user'].reshape((-1,))

        order = 4
        osr = 256
        nlev = 2
        f0 = 0.22
        Hinf = 1.25
        form = 'CRFB'

        ntf = ds.synthesizeNTF(order, osr, 2, Hinf, f0)
        a1, g1, b1, c1 = ds.realizeNTF(ntf, form)
        ABCD = ds.stuffABCD(a1, g1, b1, c1, form)

        ABCD_ref = np.array([[1., -1.6252, 0, 0, -0.0789, 0.0789],
                             [1., -0.6252, 0, 0, -0.0756, 0.0756],
                             [0, 1., 1., -1.6252, -0.2758, 0.2758],
                             [0, 1., 1., -0.6252, 0.0843, -0.0843],
                             [0, 0, 0, 1., 1., 0]])
        self.assertTrue(np.allclose(ABCD, ABCD_ref, atol=9e-5, rtol=1e-4))

        # bonus test, mapABCD - realizeNTF - stuffABCD
        a2, g2, b2, c2 = ds.mapABCD(ABCD, form)
        self.assertTrue(np.allclose(a1, a2, atol=1e-5, rtol=1e-5))
        self.assertTrue(np.allclose(g1, g2, atol=1e-5, rtol=1e-5))
        self.assertTrue(np.allclose(b1, b2, atol=1e-5, rtol=1e-5))
        self.assertTrue(np.allclose(c1, c2, atol=1e-5, rtol=1e-5))

        # We do three tests:
        # SNR from ABCD matrix
        # SNR from NTF
        # SNR from LTI obj with user specified amplitudes
        snr, amp = ds.simulateSNR(ABCD, osr, None, f0, nlev)
        self.assertTrue(np.allclose(snr, snr_ref, atol=1, rtol=5e-2))
        self.assertTrue(np.allclose(amp, amp_ref, atol=5e-1, rtol=1e-2))
        snr2, amp2 = ds.simulateSNR(ntf, osr, None, f0, nlev)
        self.assertTrue(np.allclose(snr2, snr_ref, atol=1e-5, rtol=1e-5))
        self.assertTrue(np.allclose(amp2, amp_ref, atol=1e-5, rtol=1e-5))
        amp_user = np.linspace(-100, 0, 200)[::10]
        snr_user, amp_user = ds.simulateSNR(lti(*ntf), osr=osr, amp=amp_user,
                                            f0=f0, nlev=nlev)
        self.assertTrue(np.allclose(snr_user, snr_user_ref[::10], atol=1e-5,
                                    rtol=1e-5))
        self.assertTrue(np.allclose(amp_user, amp_user_ref[::10], atol=1e-5,
                                    rtol=1e-5))

    def test_simulateSNR_2(self):
        """Test function for simulateSNR() 2/3"""
        # next test: f0 = 0
        # Load test references
        fname = pkg_resources.resource_filename(__name__,
                                                "test_data/test_snr_amp2.mat")
        amp_ref = scipy.io.loadmat(fname)['amp'].reshape((-1,))
        snr_ref = scipy.io.loadmat(fname)['snr'].reshape((-1,))
        ABCD_ref = scipy.io.loadmat(fname)['ABCD'].reshape((4, 5))

        order = 3
        osr = 256
        nlev = 2
        f0 = 0.
        Hinf = 1.25
        form = 'CIFB'

        ntf = ds.synthesizeNTF(order, osr, 2, Hinf, f0)

        a1, g1, b1, c1 = ds.realizeNTF(ntf, form)
        a1_ref = [0.008863535715733, 0.093216950269955, 0.444473912607388]
        g1_ref = [9.035620546615189e-05]
        b1_ref = [0.008863535715733, 0.093216950269955, 0.444473912607388, 1.]
        c1_ref = [1., 1., 1.]
        self.assertTrue(np.allclose(a1, a1_ref, atol=1e-9, rtol=5e-5))
        self.assertTrue(np.allclose(g1, g1_ref, atol=1e-9, rtol=5e-5))
        self.assertTrue(np.allclose(b1, b1_ref, atol=1e-9, rtol=1e-4))
        self.assertTrue(np.allclose(c1, c1_ref, atol=1e-9, rtol=2e-5))

        ABCD = ds.stuffABCD(a1, g1, b1, c1, form)
        self.assertTrue(np.allclose(ABCD, ABCD_ref, atol=9e-5, rtol=1e-4))
        snr, amp = ds.simulateSNR(ABCD, osr, None, f0, nlev)
        self.assertTrue(np.allclose(snr, snr_ref, atol=1e-5, rtol=1e-5))
        self.assertTrue(np.allclose(amp, amp_ref, atol=1e-5, rtol=1e-5))

    def test_simulateSNR_3(self):
        """Test function for simulateSNR() 3/3"""
        # next test: amp is a scalar
        fname = pkg_resources.resource_filename(__name__,
                                                "test_data/test_snr_amp2.mat")
        amp_ref = scipy.io.loadmat(fname)['amp'].reshape((-1,))[0]
        snr_ref = scipy.io.loadmat(fname)['snr'].reshape((-1,))[0]
        ABCD = scipy.io.loadmat(fname)['ABCD'].reshape((4, 5))

        order = 3
        osr = 256
        nlev = 2
        f0 = 0.
        Hinf = 1.25
        form = 'CIFB'
        ntf = ds.synthesizeNTF(order, osr, 2, Hinf, f0)

        snr, amp = ds.simulateSNR(ABCD, osr, amp_ref, f0, nlev)
        self.assertTrue(np.allclose(snr, snr_ref, atol=1e-5, rtol=1e-5))
        self.assertTrue(np.allclose(amp, amp_ref, atol=1e-5, rtol=1e-5))

