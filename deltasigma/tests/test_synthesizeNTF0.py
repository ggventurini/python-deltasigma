# -*- coding: utf-8 -*-
# test_synthesizeNTF0.py
# This module provides the tests for the synthesizeNTF0 function.
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

"""This module provides the test class for the synthesizeNTF0() function.
"""

import unittest
import numpy as np

from deltasigma._synthesizeNTF0 import synthesizeNTF0
from deltasigma._utils import cplxpair

class TestSynthesizeNTF0(unittest.TestCase):
    """Test class for synthesizeNTF0()"""

    def setUp(self):
        pass

    def test_synthesizeNTF0_1(self):
        """Test function for synthesizeNTF0() 1/8"""
        z, p, k = synthesizeNTF0(order=3, osr=64, opt=0, H_inf=1.5, f0=0.0)
        zref = [1., 1., 1.]
        pref = [.6694, .7654 + .2793j, .7654 - .2793j]
        kref = 1.
        self.assertTrue(np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF0_2(self):
        """Test function for synthesizeNTF0() 2/8"""
        # Up next: even order bandpass test
        z, p, k = synthesizeNTF0(order=4, osr=32, opt=0, H_inf=1.3, f0=.33)
        zref = [-0.4818 + 0.8763j, -0.4818 - 0.8763j, -0.4818 + 0.8763j,
                -0.4818 - 0.8763j]
        pref = [-0.5125 - 0.7018j, -0.5125 + 0.7018j, -0.3233 - 0.8240j,
                -0.3233 + 0.8240j]
        kref = 1.
        self.assertTrue(np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF0_3(self):
        """Test function for synthesizeNTF0() 3/8"""
        # repeat with zeros optimization
        z, p, k = synthesizeNTF0(order=3, osr=64, opt=1, H_inf=1.5, f0=0.0)
        zref = [1.0000 + 0.0000j, 0.9993 + 0.0380j, 0.9993 - 0.0380j]
        pref = [0.7652 - 0.2795j, 0.7652 + 0.2795j, 0.6692 + 0.0000j]
        kref = 1.
        self.assertTrue(np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF0_4(self):
        """Test function for synthesizeNTF0() 4/8"""
        # bandpass test
        z, p, k = synthesizeNTF0(order=4, osr=32, opt=1, H_inf=1.3, f0=.33)
        zref = [-0.4567 + 0.8896j, -0.4567 - 0.8896j, -0.5064 + 0.8623j,
                -0.5064 - 0.8623j]
        pref = [-0.5125 - 0.7014j, -0.5125 + 0.7014j, -0.3230 - 0.8239j,
                -0.3230 + 0.8239j]
        kref = 1.
        self.assertTrue(np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF0_5(self):
        """Test function for synthesizeNTF0() 5/8"""
        z, p, k = synthesizeNTF0(order=3, osr=64, opt=2, H_inf=1.5, f0=0.0)
        zref = [1.0000 + 0.0000j, 0.9993 + 0.0380j, 0.9993 - 0.0380j]
        pref = [0.7652 - 0.2795j, 0.7652 + 0.2795j, 0.6692 + 0.0000j]
        kref = 1.
        self.assertTrue(np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF0_6(self):
        """Test function for synthesizeNTF0() 6/8"""
        z, p, k = synthesizeNTF0(order=4, osr=32, opt=2, H_inf=1.3, f0=.33)
        zref = [-0.4818 + 0.8763j, -0.4818 - 0.8763j, -0.4818 + 0.8763j,
                -0.4818 - 0.8763j]
        pref = [-0.5125 - 0.7018j, -0.5125 + 0.7018j, -0.3233 - 0.8240j,
                -0.3233 + 0.8240j]
        kref = 1.
        self.assertTrue(np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF0_7(self):
        """Test function for synthesizeNTF0() 7/8"""
        # zeros passed explicitly
        opt = [1.0000 + 0.0000j, 0.9986 + 0.06j, 0.9986 - 0.06j,
                0.9960 + 0.0892j, 0.9960 - 0.0892j]
        z, p, k = synthesizeNTF0(order=5, osr=32, opt=opt, H_inf=1.3, f0=0.0)
        zref = [1.0000 + 0.0000j, 0.9986 + 0.06j, 0.9986 - 0.06j,
                0.9960 + 0.0892j, 0.9960 - 0.0892j]
        pref = [0.8718 - 0.0840j, 0.8718 + 0.0840j, 0.9390 - 0.1475j,
                0.9390 + 0.1475j, 0.8491 + 0.0000j]
        kref = 1.
        self.assertTrue(np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4,
                                    rtol=1e-3))
        self.assertTrue(np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4,
                                    rtol=1e-3))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeNTF0_8(self):
        """Test function for synthesizeNTF0() 8/8"""
        opt = 1
        z, p, k = synthesizeNTF0(order=4, osr=32, opt=opt, H_inf=1.3, f0=0.2)
        pref = [0.3742 - 0.7875j, 0.3742 + 0.7875j, 0.1616 - 0.8666j,
                0.1616 + 0.8666j]
        zref = [0.3358 - 0.9419j, 0.3358 + 0.9419j, 0.2819 - 0.9594j,
                0.2819 + 0.9594j]
        kref = 1.
        self.assertTrue(np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4,
                                    rtol=1e-3))
        self.assertTrue(np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4,
                                    rtol=1e-3))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

