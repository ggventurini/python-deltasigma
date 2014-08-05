# -*- coding: utf-8 -*-
# test_pulse.py
# This module provides the tests for the pulse function.
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

"""This module provides the test class for the pulse() function.
"""

import unittest
import pkg_resources
import numpy as np
import scipy.io

from nose.tools import raises

from deltasigma import pulse, axisLabels

class TestPulse(unittest.TestCase):
    """Test class for pulse()"""

    def setUp(self):
        self.H = ([1], [1, 2, 10])
        self.H0 = ([1], [1, 2, 10])
        self.H1 = ([2], [1, 2, 10])
        self.H2 = ([3], [1, 2, 10])
        fname = pkg_resources.resource_filename(__name__, "test_data/test_pulse.mat")
        self.pp2 = scipy.io.loadmat(fname)['pp']

    def test_pulse1(self):
        """Test function for pulse(): SISO 1/6"""
        # SISO
        pp = pulse([self.H], tp=(0., 1.), dt=.1, tfinal=10., nosum=False)
        pp = pp.reshape((pp.shape[0], 1))
        self.assertTrue(np.allclose(pp, self.pp2, atol=1e-6, rtol=1e-4))

    def test_pulse2(self):
        """Test function for pulse(): SIMO 2/6"""
        # SIMO
        # 1 input, 3 outputs
        pp = pulse([[self.H0, self.H1, self.H2]], tp=(0., 1.), dt=.1, tfinal=10., nosum=True)
        self.assertTrue(np.allclose(pp[:, 0, :], self.pp2, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(pp[:, 1, :], 2*self.pp2, atol=1e-6, rtol=1e-4))
        self.assertTrue(np.allclose(pp[:, 2, :], 3*self.pp2, atol=1e-6, rtol=1e-4))

    def test_pulse3(self):
        """Test function for pulse(): MISO 3/6"""
        # MISO
        # 3 inputs, one output
        # we rely here on _pulse() zipping the input tfs list
        # [H0, H1, H2] becomes [[H0], [H1], [H2]]
        pp = pulse([self.H0, self.H1, self.H2], tp=[(0., 1.)]*3, dt=.1, tfinal=10., nosum=True)
        self.assertTrue(np.allclose(pp[:, :, 0], self.pp2, atol=1e-6, rtol=1e-3))
        self.assertTrue(np.allclose(pp[:, :, 1], 2*self.pp2, atol=1e-6, rtol=1e-3))
        self.assertTrue(np.allclose(pp[:, :, 2], 3*self.pp2, atol=1e-6, rtol=1e-3))

    # FIXME ALSO CHECK MIMO TFS

    @raises(ValueError)
    def test_pulse4_tp_scalar(self):
        """Test function for pulse(): fail if tp scalar 4/6"""
        H = ((1,), (0, 0))
        pulse(H, tp=1, dt=.1, tfinal=10)

    @raises(ValueError)
    def test_pulse5_tp_not_nx2(self):
        """Test function for pulse(): fail if tp.shape not (n, 2) 5/6"""
        H = ((1,), (0, 0))
        tp = ((1, 2, 3), (4, 5, 6))
        pulse(H, tp=1, dt=.1, tfinal=10)

    @raises(ValueError)
    def test_pulse6_S_tp_mismatch(self):
        """Test function for pulse(): fail if mismatched ndac, ni 6/6"""
        H = ((1,), (0, 0))
        S = [[H], [H], [H]]
        tp = ((1, 2), (3, 4), (5, 6))
        pulse(H, tp=1, dt=.1, tfinal=10)

