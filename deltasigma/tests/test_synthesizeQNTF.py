# -*- coding: utf-8 -*-
# test_synthesizeQNTF.py
# This module provides the tests for the synthesizeQNTF function.
# Copyright 2015 Giuseppe Venturini
# This file is distributed with python-deltasigma.
#
# python-deltasigma is a 1:1 Python port of Richard Schreier's
# MATLAB delta sigma toolbox (aka "delsigma"), upon which it is heavily based.
# The delta sigma toolbox is (c) 2009, Richard Schreier.
#
# python-deltasigma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# LICENSE file for the licensing terms.
#

"""This module provides the test class for the synthesizeQNTF() function.
"""

from __future__ import division, print_function

import unittest

import numpy as np

from deltasigma import pretty_lti, synthesizeQNTF

class TestSynthesizeQNTF(unittest.TestCase):
    """Test class for synthesizeQNTF()"""

    def setUp(self):
        pass

    def test_synthesizeQNTF_1(self):
        """Test function for synthesizeQNTF() 1/4"""
        ntf = synthesizeQNTF(4, 32, 1/16, -50, -10)
        z, p, k = ntf
        p_ref = [0.573877782470855 + 0.569921695571522j,
                 0.808788367241398 + 0.002797375873482j,
                 0.591250299914031 + 0.244903892981552j,
                 0.673072277003855 - 0.278795665592338j]
        z_ref = [0.888037198535288 + 0.459771610712905j,
                 0.953044748762363 + 0.302829501298050j,
                 0.923879532511340 + 0.382683432365112j,
                 0.923879532511287 - 0.382683432365090j]
        k_ref = 1.
        np.allclose(p, p_ref, atol=1e-4, rtol=1e-3)
        np.allclose(z, z_ref, atol=1e-4, rtol=1e-3)
        np.allclose(k, k_ref, atol=1e-4, rtol=1e-3)

    def test_synthesizeQNTF_2(self):
        """Test function for synthesizeQNTF() 2/4"""
        ntf = synthesizeQNTF(4, 32, 1/16, -50, -10, 0)
        z, p, k = ntf
        z_ref = [0.910594282901269 + 0.413301405692653j,
                 0.936135618988396 + 0.351639165709982j,
                 0.888265620891369 + 0.459330150047295j,
                 0.952894107929043 + 0.303303180125290j]

        p_ref = [0.767590403998773 + 0.239841808290628j,
                 0.712362148895601 + 0.373174610786908j,
                 0.899297507221408 + 0.158408031000048j,
                 0.747910758574959 + 0.523887972745873j]
        k_ref = 1.
        np.allclose(p, p_ref, atol=1e-4, rtol=1e-3)
        np.allclose(z, z_ref, atol=1e-4, rtol=1e-3)
        np.allclose(k, k_ref, atol=1e-4, rtol=1e-3)

    def test_synthesizeQNTF_3(self):
        """Test function for synthesizeQNTF() 3/4"""
        order = 4
        osr = 32
        NG = -50
        ING = -10
        f0 = 1./16
        ntf0 = synthesizeQNTF(order, osr, f0, NG, ING)
        zeros_ref = np.array((0.8880 + 0.4598j,
                              0.9530 + 0.3028j,
                              0.9239 + 0.3827j,
                              0.9239 - 0.3827j))
        poles_ref = np.array((0.5739 + 0.5699j,
                              0.8088 + 0.0028j,
                              0.5913 + 0.2449j,
                              0.6731 - 0.2788j))
        # round down to the same precision as the results data
        # or np.sort_complex will find slightly different floats
        # to be different and sort wrongly, resulting in a failed test
        zeros = np.round(ntf0[0], 4)
        poles = np.round(ntf0[1], 4)
        assert np.allclose(zeros, zeros_ref, atol=1e-2, rtol=1e-2)
        assert np.allclose(poles, poles_ref, atol=1e-2, rtol=1e-2)
        assert ntf0[2] == 1.

    def test_synthesizeQNTF_4(self):
        """Test function for synthesizeQNTF() 4/4"""
        # order 2
        order = 2
        osr = 32
        NG = -50
        ING = -10
        f0 = 1./ 16
        ntf0 = synthesizeQNTF(order, osr, f0, NG, ING)
        z_ref = [0.900716472438935 + 0.434407454214544j,
                 0.944075182261086 + 0.329730268914908j]
        p_ref = [0.112333561599987 - 0.126178981177517j,
                 -0.00979019007164386 + 0.168653836396020j]
        k_ref = 1
        assert np.allclose(z_ref, ntf0[0])
        assert np.allclose(p_ref, ntf0[1])
        assert ntf0[2] == k_ref

