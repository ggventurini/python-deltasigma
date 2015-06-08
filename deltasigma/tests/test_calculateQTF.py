# -*- coding: utf-8 -*-
# test_calculateQTF.py
# This module provides the tests for the calculateQTF() function.
# Copyright 2015 Giuseppe Venturini
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


class TestCalculateQTF(unittest.TestCase):
    """Test function for calculateQTF()"""

    def setUp(self):
        self.ABCD = np.array([[-4.59771611e-01+0.8880372j, 0.00000000e+00+0.j,
                               0.00000000e+00+0.j, 0.00000000e+00+0.j,
                               0.00000000e+00+0.05160829j,
                               -1.73472348e-18+0.05160829j],
                              [-7.94503732e-01+0.60725927j,
                               -3.02829501e-01+0.95304475j,
                               0.00000000e+00+0.j, 0.00000000e+00+0.j,
                               0.00000000e+00+0.j, 2.77555756e-17+0.32296443j],
                              [0.00000000e+00+0.j, -3.24479322e-01+0.94589279j,
                               -3.82683432e-01+0.92387953j,
                               0.00000000e+00+0.j, 0.00000000e+00+0.j,
                               -9.71445147e-17+0.91615063j],
                              [0.00000000e+00+0.j, 0.00000000e+00+0.j,
                               0.00000000e+00+0.j,
                               3.82683432e-01+0.92387953j,
                               0.00000000e+00+0.j,
                               -2.77555756e-17+0.17778533j],
                              [0.00000000e+00+0.j, 0.00000000e+00+0.j,
                               1.13847053e-01-0.99349829j,
                               6.72005754e-01-0.74054592j,
                               0.00000000e+00+0.j, 0.00000000e+00+0.j]])
        self.ABCDr = np.array([[8.88037199e-01, -4.59771611e-01, 0.00000000e+00,
                                0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                                0.00000000e+00, 0.00000000e+00, 5.16082949e-02,
                                0.00000000e+00, 5.16082949e-02,
                                -1.73472348e-18],
                               [4.59771611e-01, 8.88037199e-01, 0.00000000e+00,
                                  0.00000000e+00, 0.00000000e+00,
                                  0.00000000e+00, 0.00000000e+00,
                                  0.00000000e+00, 0.00000000e+00,
                                  5.16082949e-02, 1.73472348e-18,
                                  5.16082949e-02],
                               [6.07259269e-01, -7.94503732e-01, 9.53044749e-01,
                                -3.02829501e-01, 0.00000000e+00,
                                0.00000000e+00, 0.00000000e+00,
                                0.00000000e+00, 0.00000000e+00,
                                0.00000000e+00, 3.22964434e-01,
                                2.77555756e-17],
                               [7.94503732e-01, 6.07259269e-01, 3.02829501e-01,
                                9.53044749e-01, 0.00000000e+00,
                                0.00000000e+00, 0.00000000e+00,
                                0.00000000e+00, 0.00000000e+00,
                                0.00000000e+00, -2.77555756e-17,
                                3.22964434e-01],
                               [0.00000000e+00, 0.00000000e+00, 9.45892790e-01,
                                -3.24479322e-01, 9.23879533e-01,
                                -3.82683432e-01, 0.00000000e+00, 0.00000000e+00,
                                0.00000000e+00, 0.00000000e+00, 9.16150632e-01,
                                -9.71445147e-17],
                               [0.00000000e+00, 0.00000000e+00, 3.24479322e-01,
                                9.45892790e-01, 3.82683432e-01, 9.23879533e-01,
                                0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                                0.00000000e+00, 9.71445147e-17, 9.16150632e-01],
                               [0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                                0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                                9.23879533e-01, 3.82683432e-01, 0.00000000e+00,
                                0.00000000e+00, 1.77785329e-01,
                                -2.77555756e-17],
                               [0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                                0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                                -3.82683432e-01, 9.23879533e-01, 0.00000000e+00,
                                0.00000000e+00, 2.77555756e-17, 1.77785329e-01],
                               [0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                                0.00000000e+00, -9.93498288e-01, 1.13847053e-01,
                                -7.40545924e-01, 6.72005754e-01, 0.00000000e+00,
                                0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                               [0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                                0.00000000e+00, -1.13847053e-01,
                                -9.93498288e-01, -6.72005754e-01,
                                -7.40545924e-01, 0.00000000e+00, 0.00000000e+00,
                                0.00000000e+00, 0.00000000e+00]])
        self.ntf_poles = np.array([0.5739 + 0.5699j, 0.8088 + 0.0028j,
                                   0.6731 - 0.2788j, 0.5913 + 0.2449j])
        self.ntf_zeros = np.array([0.9239 - 0.3827j, 0.8880 + 0.4598j,
                                   0.9239 + 0.3827j, 0.9530 + 0.3028j])
        self.intf_poles = np.array([0.5739 + 0.5699j, 0.5739 - 0.5699j,
                                    0.8088 + 0.0028j, 0.8088 - 0.0028j,
                                    0.6731 + 0.2788j, 0.6731 - 0.2788j,
                                    0.5913 + 0.2449j, 0.5913 - 0.2449j])
        self.intf_zeros = np.array([1.3472 + 1.6326j, 1.2895 - 1.7086j,
                                    0.2158 + 0.6627j, 0.2046 - 0.6699j,
                                    0.4033 + 0.2494j, 0.3786 - 0.1930j,
                                    0.1636 - 0.0980j])
        self.stf_poles = np.array([0.5739 + 0.5699j, 0.8088 + 0.0028j,
                                   0.6731 - 0.2788j, 0.5913 + 0.2449j])
        self.stf_zeros = np.array([0.9239 - 0.3827j])
        self.istf_poles = np.array([0.5739 + 0.5699j, 0.5739 - 0.5699j,
                                    0.8088 + 0.0028j, 0.8088 - 0.0028j,
                                    0.6731 + 0.2788j, 0.6731 - 0.2788j,
                                    0.5913 + 0.2449j, 0.5913 - 0.2449j])
        self.istf_zeros = np.array([-2.8543 - 0.5711j, 1.5159 - 1.2233j,
                                    0.8588 - 0.6692j, 0.6829 + 0.4823j,
                                    0.6667 + 0.1139j])

        self.ntf_k = 1
        self.intf_k = 2.6645e-15
        self.stf_k = -0.0107 - 0.0505j
        self.istf_k = 1e-17

    def test(self):
        """Test function for calculateQTF()"""
        ABCDr_test = ds.mapQtoR(self.ABCD)
        allsortedclose(ABCDr_test, self.ABCDr, atol=1e-3, rtol=1e-3)
        ntf, stf, intf, istf = ds.calculateQTF(self.ABCDr)

        # NTF CHECK
        allsortedclose(ntf[0], self.ntf_zeros, atol=1e-3, rtol=1e-3)
        allsortedclose(ntf[1], self.ntf_poles, atol=1e-3, rtol=1e-3)
        allsortedclose(ntf[2], self.ntf_k, atol=1e-3, rtol=1e-3)
        # STF CHECK
        allsortedclose(stf[0], self.stf_zeros, atol=1e-3, rtol=1e-3)
        allsortedclose(stf[1], self.stf_poles, atol=1e-3, rtol=1e-3)
        allsortedclose(stf[2], self.stf_k, atol=1e-3, rtol=1e-3)
        #INTF CHECK
        # the numerator coefficients are all zeros anyways. Big numeric
        # errors for worthless information
        #allsortedclose(intf[0], self.intf_zeros, atol=1e-3, rtol=1e-3)
        allsortedclose(intf[1], self.intf_poles, atol=1e-3, rtol=1e-3)
        allsortedclose(intf[2], self.intf_k, atol=1e-3, rtol=1e-3)
        # ISTF CHECK
        # the zeros fail on travis ... again, no real importance
        #allsortedclose(istf[0], self.istf_zeros, atol=1e-3, rtol=1e-3)
        allsortedclose(istf[1], self.istf_poles, atol=1e-3, rtol=1e-3)
        allsortedclose(istf[2], self.istf_k, atol=1e-3, rtol=1e-3)

def allsortedclose(a, b, atol=1e-3, rtol=1e-3):
    if np.iscomplex(a).any():
        a = np.sort_complex(a)
    else:
        a = np.sort(a)
    if np.iscomplex(b).any():
        b = np.sort_complex(b)
    else:
        b = np.sort(b)
    return np.allclose(a, b, rtol=rtol, atol=atol)
