# -*- coding: utf-8 -*-
# test_synthesizeNTF.py
# This module provides the tests for the synthesizeNTF function.
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

"""This module provides the test class for the synthesizeNTF() function.
"""

import unittest
import numpy as np

from nose.tools import raises

from deltasigma import synthesizeNTF
from deltasigma._utils import cplxpair as cpx

class TestSynthesizeNTF(unittest.TestCase):
    """Test class for synthesizeNTF()"""

    def setUp(self):
        pass

    def test_synthesizeNTF_1(self):
        """Test function for synthesizeNTF() 1/15"""
        # synthesizeNTF should have as default values:
        # order=3, osr=64, opt=0, H_inf=1.5, f0=0.0
        z, p, k = synthesizeNTF()
        zref = [1., 1., 1.]
        pref = [.6694, .7654 + .2793j, .7654 - .2793j]
        kref = 1.
        self.assertTrue(np.allclose(cpx(z), cpx(zref), atol=1e-4, rtol=1e-4))
        self.assertTrue( np.allclose(cpx(p), cpx(pref), atol=1e-4, rtol=1e-4))
        self.assertTrue( np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF_2(self):
        """Test function for synthesizeNTF() 2/15"""
        # Up next: bandpass test
        z, p, k = synthesizeNTF(order=4, osr=32, opt=0, H_inf=1.3, f0=.33)
        zref = [-0.4818 + 0.8763j, -0.4818 - 0.8763j, -0.4818 + 0.8763j,
                -0.4818 - 0.8763j]
        pref = [-0.5125 - 0.7018j, -0.5125 + 0.7018j, -0.3233 - 0.8240j, 
                -0.3233 + 0.8240j]
        kref = 1.
        self.assertTrue(np.allclose(cpx(z), cpx(zref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(cpx(p), cpx(pref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF_3(self):
        """Test function for synthesizeNTF() 3/15"""
        # repeat with zeros optimization
        z, p, k = synthesizeNTF(opt=1)
        zref = [1.0000 + 0.0000j, 0.9993 + 0.0380j, 0.9993 - 0.0380j]
        pref = [0.7652 - 0.2795j, 0.7652 + 0.2795j, 0.6692 + 0.0000j]
        kref = 1.
        self.assertTrue(np.allclose(cpx(z), cpx(zref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(cpx(p), cpx(pref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF_4(self):
        """Test function for synthesizeNTF() 4/15"""
        z, p, k = synthesizeNTF(order=4, osr=32, opt=1, H_inf=1.3, f0=.33)
        zref = [-0.4567 + 0.8896j, -0.4567 - 0.8896j, -0.5064 + 0.8623j,
                -0.5064 - 0.8623j]
        pref = [-0.5125 - 0.7014j, -0.5125 + 0.7014j, -0.3230 - 0.8239j,
                -0.3230 + 0.8239j]
        kref = 1.
        self.assertTrue(np.allclose(cpx(z), cpx(zref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(cpx(p), cpx(pref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF_5(self):
        """Test function for synthesizeNTF() 5/15"""
        z, p, k = synthesizeNTF(opt=2)
        zref = [1.0000 + 0.0000j, 0.9993 + 0.0380j, 0.9993 - 0.0380j]
        pref = [0.7652 - 0.2795j, 0.7652 + 0.2795j, 0.6692 + 0.0000j]
        kref = 1.
        self.assertTrue(np.allclose(cpx(z), cpx(zref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(cpx(p), cpx(pref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF_6(self):
        """Test function for synthesizeNTF() 6/15"""
        z, p, k = synthesizeNTF(order=4, osr=32, opt=2, H_inf=1.3, f0=.33)
        zref = [-0.4818 + 0.8763j, -0.4818 - 0.8763j, -0.4818 + 0.8763j,
                -0.4818 - 0.8763j]
        pref = [-0.5125 - 0.7018j, -0.5125 + 0.7018j, -0.3233 - 0.8240j,
                -0.3233 + 0.8240j]
        kref = 1.
        self.assertTrue(np.allclose(cpx(z), cpx(zref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(cpx(p), cpx(pref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF_7(self):
        """Test function for synthesizeNTF() 7/15"""
        # opt = 3
        z, p, k = synthesizeNTF(opt=3)
        zref = [1.0000 + 0.0000j, 0.9993 + 0.0380j, 0.9993 - 0.0380j]
        pref = [0.7652 - 0.2795j, 0.7652 + 0.2795j, 0.6692 + 0.0000j]
        kref = 1.
        self.assertTrue(np.allclose(cpx(z), cpx(zref), atol=1e-4, rtol=1e-3))
        self.assertTrue(np.allclose(cpx(p), cpx(pref), atol=1e-4, rtol=1e-3))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF_8(self):
        """Test function for synthesizeNTF() 8/15"""
        z, p, k = synthesizeNTF(order=4, osr=32, opt=3, H_inf=1.3, f0=.33)
        zref = [-0.4567 + 0.8896j, -0.4567 - 0.8896j, -0.5064 + 0.8623j,
                -0.5064 - 0.8623j]
        pref = [-0.5125 - 0.7014j, -0.5125 + 0.7014j, -0.3230 - 0.8239j,
                -0.3230 + 0.8239j]
        kref = 1.
        self.assertTrue(np.allclose(cpx(z), cpx(zref), atol=1e-4, rtol=1e-3))
        self.assertTrue(np.allclose(cpx(p), cpx(pref), atol=1e-4, rtol=1e-3))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF_9(self):
        """Test function for synthesizeNTF() 9/15"""
        # opt = 4
        z, p, k = synthesizeNTF(opt=4)
        zref = [1.0000 + 0.0000j, 0.9993 + 0.0380j, 0.9993 - 0.0380j]
        pref = [0.7652 - 0.2795j, 0.7652 + 0.2795j, 0.6692 + 0.0000j]
        kref = 1.
        self.assertTrue(np.allclose(cpx(z), cpx(zref), atol=1e-4, rtol=1e-3))
        self.assertTrue(np.allclose(cpx(p), cpx(pref), atol=1e-4, rtol=1e-3))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF_10(self):
        """Test function for synthesizeNTF() 10/15"""
        # Up next: odd order lowpass w zero at center band test
        z, p, k = synthesizeNTF(order=5, osr=32, opt=4, H_inf=1.3, f0=0.0)
        zref = [1.0000 + 0.0000j, 0.9986 + 0.0531j, 0.9986 - 0.0531j,
                0.9960 + 0.0892j, 0.9960 - 0.0892j]
        pref = [0.8718 - 0.0840j, 0.8718 + 0.0840j, 0.9390 - 0.1475j,
                0.9390 + 0.1475j, 0.8491 + 0.0000j]
        kref = 1.
        self.assertTrue(np.allclose(cpx(z), cpx(zref), atol=1e-4, rtol=1e-3))
        self.assertTrue(np.allclose(cpx(p), cpx(pref), atol=1e-4, rtol=1e-3))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF_11(self):
        """Test function for synthesizeNTF() 11/15"""
        # zeros passed explicitly
        opt = [1.0000 + 0.0000j, 0.9986 + 0.06j, 0.9986 - 0.06j,
                0.9960 + 0.0892j, 0.9960 - 0.0892j]
        z, p, k = synthesizeNTF(order=5, osr=32, opt=opt, H_inf=1.3, f0=0.0)
        zref = [1.0000 + 0.0000j, 0.9986 + 0.06j, 0.9986 - 0.06j,
                0.9960 + 0.0892j, 0.9960 - 0.0892j]
        pref = [0.8718 - 0.0840j, 0.8718 + 0.0840j, 0.9390 - 0.1475j,
                0.9390 + 0.1475j, 0.8491 + 0.0000j]
        kref = 1.
        self.assertTrue(np.allclose(cpx(z), cpx(zref), atol=1e-4, rtol=1e-3))
        self.assertTrue(np.allclose(cpx(p), cpx(pref), atol=1e-4, rtol=1e-3))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    @raises(ValueError)
    def test_synthesizeNTF_12(self):
        """Test function for synthesizeNTF() 12/15"""
        # f0 > 0.5 -> ValueError
        z, p, k = synthesizeNTF(order=2, osr=32, opt=0, H_inf=1.3, f0=0.7)

    @raises(ValueError)
    def test_synthesizeNTF_13(self):
        """Test function for synthesizeNTF() 13/15"""
        # f0 > 0. and order odd -> ValueError
        z, p, k = synthesizeNTF(order=3, osr=32, opt=0, H_inf=1.3, f0=0.3)

    @raises(ValueError)
    def test_synthesizeNTF_14(self):
        """Test function for synthesizeNTF() 14/15"""
        # 1 < len(opt) < order
        z, p, k = synthesizeNTF(order=3, osr=32, opt=[0., 0.], H_inf=1.3, f0=0.3)

    @raises(ValueError)
    def test_synthesizeNTF_15(self):
        """Test function for synthesizeNTF() 15/15"""
        # order < len(opt)
        z, p, k = synthesizeNTF(order=3, osr=32, opt=[0., 0., 0., 0.], H_inf=1.3, f0=0.)

