# -*- coding: utf-8 -*-
# test_synthesizeChebyshevNTF.py
# This module provides the tests for the synthesizeChebyshevNTF function.
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

"""This module provides the test class for the synthesizeChebyshevNTF() function.
"""

import unittest
import numpy as np
import deltasigma as ds

from warnings import catch_warnings

from nose.tools import raises

from deltasigma._utils import cplxpair

class TestSynthesizeChebyshevNTF(unittest.TestCase):
    """Test class for synthesizeChebyshevNTF()"""

    def setUp(self):
        pass

    def test_synthesizeChebyshevNTF_1(self):
        """Test function for synthesizeChebyshevNTF() 1/3"""
        z, p, k = ds.synthesizeChebyshevNTF()
        zref = [1., .9991 + 0.0425j, .9991 - 0.0425j]
        pref = [.6609, .7686 + .2858j, .7686 - .2858j]
        kref = 1.
        self.assertTrue(np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    def test_synthesizeChebyshevNTF_2(self):
        """Test function for synthesizeChebyshevNTF() 2/3"""
        with catch_warnings(record=True) as w:
            z, p, k = ds.synthesizeChebyshevNTF(order=4, OSR=32, opt=1, H_inf=1.5, f0=.33)
            self.assertTrue(len(w) > 0)
        zref = [-.4513 + .8924j, -.4513 - .8924j, -.5122 + 0.8589j, -.5122 - 0.8589j]
        pref = [-.2249 + .7665j, -.2249 - .7665j, -.5506 + .6314j, -.5506 - .6314j]
        kref = 1.
        self.assertTrue(np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-4))
        self.assertTrue(np.allclose(k, kref, atol=1e-4, rtol=1e-4))

    @raises(ValueError)
    def test_synthesizeChebyshevNTF_3(self):
        """Test function for synthesizeChebyshevNTF() 3/3"""
        z, p, k = ds.synthesizeChebyshevNTF(order=5, OSR=32, opt=0, H_inf=1.5, f0=.33)

