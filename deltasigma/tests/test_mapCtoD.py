# -*- coding: utf-8 -*-
# test_mapCtoD.py
# This module provides the tests for the mapCtoD function.
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

"""This module provides the test class for the mapCtoD() function.
"""

from __future__ import division, print_function

import unittest
import numpy as np
import deltasigma as ds

from nose.tools import raises

class TestMapCtoD(unittest.TestCase):
    """Test class for mapCtoD()"""

    def setUp(self):
        pass

    def test_mapCtoD_1(self):
        """Test function for mapCtoD() 1/6"""
        # The first test comes straight from the DSToolbox doc
        # We map the standard continuous time DS modulator to its DT counterpart
        # and we check whether the TF is (1 - 1/z)**2
        LFc = (np.array([[0., 0.], [1., 0.]]),
               np.array([[1., -1.], [0., -1.5]]),
               np.array([[0., 1.]]),
               np.array(([[0., 0.]]))
              )
        tdac = [0, 1]
        LF, Gp = ds.mapCtoD(LFc, tdac)
        ABCD = np.vstack((np.hstack((LF[0], LF[1])),
                          np.hstack((LF[2], LF[3]))
                        ))
        H = ds.calculateTF(ABCD)
        self.assertTrue(np.allclose(H[0].num, [ 1., -2.,  1.], atol=1e-8, rtol=1e-5))
        self.assertTrue(np.allclose(H[0].den, [1., 0., 0.], atol=1e-8, rtol=1e-5))

    def test_mapCtoD_2(self):
        """Test function for mapCtoD() 2/6"""
        # zero delay NRZ test
        ABCDc = np.array([[0., 0.,  0., 0.04440879, -0.04440879],
                          [1., 0., -0.00578297, 0., -0.23997611],
                          [0., 1.,  0., 0., -0.67004646],
                          [0., 0.,  1., 0.,  0.]])
        tdac = np.array([[-1., -1], [0., 1.]])
        sys_d, Gp = ds.mapCtoD(ABCDc, tdac)
        ABCD = np.vstack((np.hstack((sys_d[0], sys_d[1])),
                          np.hstack((sys_d[2], sys_d[3]))
                        ))
        ABCDref = np.array([[1., 0., 0., 0.04440879, -0.04440879],
                            [0.99903645, 0.99710991, -0.0057774, 0.0221937, -0.26000208],
                            [0.49975909, 0.99903645,  0.99710991, 0.00739932, -0.79673041],
                            [0., 0., 1., 0., 0.]])
        self.assertTrue(np.allclose(ABCD, ABCDref, atol=1e-8, rtol=1e-5))

    def test_mapCtoD_3(self):
        """Test function for mapCtoD() 3/6"""
        # .1 delay NRZ
        ABCDc = np.array([[0., 0., 0., 0.0444, -0.0444],
                          [1., 0., -0.0058, 0., -0.2440],
                          [0., 1.,  0., 0., -0.6942],
                          [0., 0.,  1., 0., -0.0682]])

        tdac = np.array([[-1., -1], [0.1, 1.1]])
        sys_d, Gp = ds.mapCtoD(ABCDc, tdac)
        ABCD = np.vstack((np.hstack((sys_d[0], sys_d[1])),
                          np.hstack((sys_d[2], sys_d[3]))
                        ))
        ABCDref = np.array([[1., 0., 0., -0.0044, 0.0444, -0.0400],
                            [0.999, 0.9971, -0.0058, -0.0282, 0.0222, -0.2358],
                            [0.4998, 0.999, 0.9971, -0.0944, 0.0074, -0.7285],
                            [0., 0., 0., 0., 0., 1.],
                            [0., 0., 1., -0.0682, 0., 0.]])
        self.assertTrue(np.allclose(ABCD, ABCDref, atol=1e-4, rtol=1e-4))

    def test_mapCtoD_4(self):
        """Test function for mapCtoD() 4/6"""
        # Non-zero f0 test
        f0 = .2
        ABCDc = np.array([[0., 0., 1., -1.],
                          [1., 0., 0., -1.5],
                          [0., 1., 0., 0.]])
        tdac = np.array([0.1, 1.1])
        sys_d, Gp = ds.mapCtoD(ABCDc, tdac, f0=f0)
        ABCD = np.vstack((np.hstack((sys_d[0], sys_d[1])),
                          np.hstack((sys_d[2], sys_d[3]))
                        ))
        ABCDref = np.array([[1., 0., -0.1, 0.9355, -0.9],
                            [1., 1., -0.245, 0.4784, -1.7550],
                            [0., 0., 0., 0., 1.],
                            [0., 1., 0., 0., 0]])
        self.assertTrue(np.allclose(ABCD, ABCDref, atol=1e-4, rtol=1e-4))

    def test_mapCtoD_5(self):
        """Test function for mapCtoD() 5/6"""
        # Non-zero f0 test, short DAC pulse
        f0 = .2
        ABCDc = np.array([[0., 0., 1., -1.],
                          [1., 0., 0., -1.5],
                          [0, 1, 0, 0]])
        tdac = np.array([0.1, 0.4])
        sys_d, Gp = ds.mapCtoD(ABCDc, tdac, f0=f0)
        ABCD = np.vstack((np.hstack((sys_d[0], sys_d[1])),
                          np.hstack((sys_d[2], sys_d[3]))
                        ))
        print(ABCD)
        ABCDref = np.array([[1., 0, 0.9355, -0.3000],
                            [1., 1., 0.4784, -0.6750],
                            [0., 1., 0., 0.]])
        self.assertTrue(np.allclose(ABCD, ABCDref, atol=1e-4, rtol=1e-4))

    @raises(ValueError)
    def test_mapCtoD_6(self):
        """Test function for mapCtoD() 6/6"""
        # wrong ni
        f0 = 0.0
        ABCDc = np.array([[0., 0., 1., -1.],
                          [1., 0., 0., -1.5],
                          [0, 1, 0, 0]])
        tdac = np.array([0.1, 0.4, 0.1, 3.])
        sys_d, Gp = ds.mapCtoD(ABCDc, tdac, f0=f0)

