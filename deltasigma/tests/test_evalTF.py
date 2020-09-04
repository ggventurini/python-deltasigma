# -*- coding: utf-8 -*-
# test_evalTF.py
# This module provides the tests for the evalTF function.
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

from scipy.signal import lti, tf2zpk
from deltasigma._utils import empty

from nose.tools import raises

class TestEvalTF(unittest.TestCase):
    """Test functions for evalTF()"""
    def setUp(self):
        num = np.poly([3, 0.3, 1])
        den = np.poly([2, 0.5, .25])
        H = (num, den)
        tstr1 = empty()
        (tstr1.form, tstr1.num, tstr1.den) = ('coeff', num, den)
        tstr2 = empty()
        tstr2.form = 'zp'
        (tstr2.zeros, tstr2.poles, tstr2.gain) = tf2zpk(num, den)
        z = np.exp(1j * np.linspace(0, 2*np.pi, num=129, endpoint=True))
        self.h1 = ds.evalTF(tstr1, z)
        self.h2 = ds.evalTF(tstr2, z)
        self.h3 = ds.evalTF(H, z)
        self.h4 = ds.evalTF(lti(tstr2.zeros, tstr2.poles, tstr2.gain), z)
        h5tf = lti(tstr2.zeros, tstr2.poles, tstr2.gain).to_ss()
        self.h5 = ds.evalTF((h5tf.A, h5tf.B, h5tf.C, h5tf.D), z)
        h6tf = np.vstack((np.hstack((h5tf.A, h5tf.B)),
                          np.hstack((h5tf.C, np.atleast_2d(h5tf.D)))))
        self.h6 = ds.evalTF(h6tf, z)

    def test_evalTF_first(self):
        """Test evalTF num/den obj form"""
        self.assertTrue(np.allclose(
            np.abs(self.h1), np.abs(self.h2), atol=1e-8, rtol=1e-5))

    def test_evalTF_second(self):
        """Test evalTF zp obj form"""
        self.assertTrue(np.allclose(self.h1, self.h3, atol=1e-8, rtol=1e-5))

    def test_evalTF_third(self):
        """Test evalTF zpk and num/den tuple form"""
        self.assertTrue(np.allclose(self.h3, self.h4, atol=1e-8, rtol=1e-5))

    def test_evalTF_fourth(self):
        """Test evalTF ABCD tuple form"""
        self.assertTrue(np.allclose(self.h5, self.h4, atol=1e-8, rtol=1e-5))

    def test_evalTF_fifth(self):
        """Test evalTF ABCD matrix form"""
        self.assertTrue(np.allclose(self.h6, self.h5, atol=1e-8, rtol=1e-5))

    @raises(TypeError)
    def test_evalTF_sixth(self):
        """Test unknown forms error"""
        ds.evalTF(empty(), [1, 2, 3])
