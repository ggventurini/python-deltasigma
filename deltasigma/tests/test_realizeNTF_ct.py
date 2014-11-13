# -*- coding: utf-8 -*-
# test_realizeNTF_ct.py
# This module provides the tests for the realizeNTF_ct function.
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

"""This module provides the test class for the realizeNTF_ct() function.
"""

import unittest
import numpy as np
import deltasigma as ds

from nose.tools import raises

class TestRealizeNTF_CT(unittest.TestCase):
    """Test class for realizeNTF_ct()"""

    def setUp(self):
        # FB NRTZ 0 delay tdac default
        self.ntf = (np.array([1., 1.]), np.array([0., 0.]), 1)
        self.ABCDc_ref1 = np.array(((0, 0, 1, -1),
                                   (1, 0, 0, -1.5),
                                   (0, 1, 0, 0)), dtype=np.float)
        self.tdac2_ref1 = np.array(((-1, -1),
                                    (0,  1)))
        # FB NRTZ 0 delay tdac in short list form
        self.ABCDc_ref2 = np.array(((0, 0, 1, -1),
                                    (1, 0, 0, -1.5),
                                    (0, 1, 0, 0)), dtype=np.float)
        self.tdac2_ref2 = np.array(((-1, -1),
                                    (0,  1)))
        # FB multi-timing DAC
        self.ntf3 = (np.array([1., 1., 1]), np.array([0., 0., 0]), 1)
        self.tdac3 = [[1, 2], [1, 2], [[0.5, 1], [1, 1.5]], []]
        self.ABCDc_ref3 = np.array(((0., 0., 0., 1., -1., 0., 0., 0.),
                                    (1., 0., 0., 0., 0., -3., 0., 0.),
                                    (0., 1., 0., 0., 0., 0., -6., -2.6667),
                                    (0., 0., 1., 0., 0., 0., 0., 0.)))
        self.tdac2_ref3 = np.array(((-1., -1.),
                                    (1., 2.),
                                    (1., 2.),
                                    (0.5, 1.),
                                    (1., 1.5)))
        # test for FF
        self.ABCDc_ref4 = np.array([[0., 0., 1., -1.],
                                    [1., 0., 0.,  0.],
                                    [ 1.5, 1., 0., 0.]])
        self.tdac2_ref4 = np.array([[-1., -1.], [ 0., 1.]])
        # test for FF - NRZ DAC pulse delayed by .5
        self.ABCDc_ref5 = np.array([[0., 0., 1., -1.],
                                    [1., 0., 0., 0.],
                                    [2., 1., 0., -0.8750]])
        self.tdac2_ref5 = np.array([[-1., -1.], [ 0.5, 1.5]])

    def test_realizeNTF_ct_1(self):
        """Test function for realizeNTF_ct() 1/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf, 'FB')
        self.assertTrue(np.allclose(ABCDc, self.ABCDc_ref1, atol=1e-8,
                                    rtol=1e-5))
        self.assertTrue(np.allclose(tdac2, self.tdac2_ref1, atol=1e-8,
                                    rtol=1e-5))

    def test_realizeNTF_ct_2(self):
        """Test function for realizeNTF_ct() 2/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf, 'FB', tdac=[0, 1])
        self.assertTrue(np.allclose(ABCDc, self.ABCDc_ref2, atol=1e-8,
                                    rtol=1e-5))
        self.assertTrue(np.allclose(tdac2, self.tdac2_ref2, atol=1e-8,
                                    rtol=1e-5))

    def test_realizeNTF_ct_3(self):
        """Test function for realizeNTF_ct() 3/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf3, 'FB', self.tdac3)
        self.assertTrue(np.allclose(ABCDc, self.ABCDc_ref3, atol=1e-8,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(tdac2, self.tdac2_ref3, atol=1e-8,
                                    rtol=1e-4))

    def test_realizeNTF_ct_4(self):
        """Test function for realizeNTF_ct() 4/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf, 'FF')
        self.assertTrue(np.allclose(ABCDc, self.ABCDc_ref4, atol=1e-8,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(tdac2, self.tdac2_ref4, atol=1e-8,
                                    rtol=1e-4))

    def test_realizeNTF_ct_5(self):
        """Test function for realizeNTF_ct() 5/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf, 'FF', tdac=[.5, 1.5])
        self.assertTrue(np.allclose(ABCDc, self.ABCDc_ref5, atol=1e-8,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(tdac2, self.tdac2_ref5, atol=1e-8,
                                    rtol=1e-4))

    @raises(ValueError)
    def test_realizeNTF_ct_6(self):
        """Test function for realizeNTF_ct() 6/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf, 'DUMMY', tdac=[.5, 1.5])

    @raises(ValueError)
    def test_realizeNTF_ct_7(self):
        """Test function for realizeNTF_ct() 7/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf, 'FF', tdac=[.5, 1.5, 3.])

    @raises(ValueError)
    def test_realizeNTF_ct_8(self):
        """Test function for realizeNTF_ct() 8/15"""
        tdac= [[-1, -1], [.5, 1.5], [.5, .8], [0.1, 2.]]
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf3, 'FB', tdac)

    @raises(ValueError)
    def test_realizeNTF_ct_9(self):
        """Test function for realizeNTF_ct() 9/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf3, 'FF', self.tdac3)

    def test_realizeNTF_ct_10(self):
        """Test function for realizeNTF_ct() 10/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf, 'FB', method='NTF')
        self.assertTrue(np.allclose(ABCDc, self.ABCDc_ref1, atol=1e-8,
                                    rtol=1e-5))
        self.assertTrue(np.allclose(tdac2, self.tdac2_ref1, atol=1e-8,
                                    rtol=1e-5))

    def test_realizeNTF_ct_11(self):
        """Test function for realizeNTF_ct() 11/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf, 'FB', tdac=[0, 1], method='NTF')
        self.assertTrue(np.allclose(ABCDc, self.ABCDc_ref2, atol=1e-8,
                                    rtol=1e-5))
        self.assertTrue(np.allclose(tdac2, self.tdac2_ref2, atol=1e-8,
                                    rtol=1e-5))

    def test_realizeNTF_ct_12(self):
        """Test function for realizeNTF_ct() 12/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf3, 'FB', self.tdac3, method='NTF')
        self.assertTrue(np.allclose(ABCDc, self.ABCDc_ref3, atol=1e-8,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(tdac2, self.tdac2_ref3, atol=1e-8,
                                    rtol=1e-4))

    def test_realizeNTF_ct_13(self):
        """Test function for realizeNTF_ct() 13/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf, 'FF', method='NTF')
        self.assertTrue(np.allclose(ABCDc, self.ABCDc_ref4, atol=1e-8,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(tdac2, self.tdac2_ref4, atol=1e-8,
                                    rtol=1e-4))

    def test_realizeNTF_ct_14(self):
        """Test function for realizeNTF_ct() 14/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf, 'FF', tdac=[.5, 1.5],
                                        method='NTF')
        self.assertTrue(np.allclose(ABCDc, self.ABCDc_ref5, atol=1e-8,
                                    rtol=1e-4))
        self.assertTrue(np.allclose(tdac2, self.tdac2_ref5, atol=1e-8,
                                    rtol=1e-4))

    @raises(ValueError)
    def test_realizeNTF_ct_15(self):
        """Test function for realizeNTF_ct() 15/15"""
        ABCDc, tdac2 = ds.realizeNTF_ct(self.ntf3, 'FB', method='DUMMY')

