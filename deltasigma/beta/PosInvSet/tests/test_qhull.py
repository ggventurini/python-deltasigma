# -*- coding: utf-8 -*-
# test_qhull.py
# This module provides the tests for the qhull() function.
# Copyright 2020 Yuki Fukuda
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

from nose.tools import raises

class TestQHull(unittest.TestCase):
    """Test function for qhull()"""

    def setUp(self):
        self.points = np.array([[0.20184693, 0.97717224],
                             [0.20297392, 0.4782398 ],
                             [0.98686318, 0.68111812],
                             [0.10370579, 0.68252684],
                             [0.24010152, 0.06761969],
                             [0.37566978, 0.26573704],
                             [0.99609644, 0.45377267],
                             [0.56771729, 0.23667973],
                             [0.80680514, 0.02321781],
                             [0.44561283, 0.78933871],
                             [0.79447218, 0.70143051],
                             [0.99636547, 0.60034883],
                             [0.50759036, 0.65499842],
                             [0.47678139, 0.72627565],
                             [0.52188126, 0.61206678],
                             [0.17915326, 0.19261078],
                             [0.1540633 , 0.84502877],
                             [0.58649626, 0.19831919],
                             [0.02067492, 0.90157589],
                             [0.62496933, 0.58129522],
                             [0.48413661, 0.11267991],
                             [0.54166052, 0.70674753],
                             [0.82643576, 0.57895598],
                             [0.80705596, 0.41159023],
                             [0.89475419, 0.85872548],
                             [0.23566104, 0.17886997],
                             [0.97317125, 0.24226311],
                             [0.47825798, 0.91981332],
                             [0.00408167, 0.71863827],
                             [0.93084148, 0.65450774]])

    def test_qhull_generation(self):
        V, E, N, O = ds.qhull(self.points)
        vert = np.array([ 4,  8, 26,  6, 11,  2, 24,  0, 18, 28, 15], dtype='int32')
        edge = np.array([[ 4,  8],
                         [26,  8],
                         [15, 28],
                         [15,  4],
                         [18, 28],
                         [18,  0],
                         [24,  0],
                         [ 6, 11],
                         [ 6, 26],
                         [ 2, 11],
                         [ 2, 24]], dtype='int32')
        normals = np.array([[-0.07811175, -0.99694461],
                            [ 0.79635185, -0.60483364],
                            [-0.94882979, -0.31578794],
                            [-0.89883361, -0.43829002],
                            [-0.99591157,  0.09033356],
                            [-0.38508417,  0.92288146],
                            [ 0.16849761,  0.98570206],
                            [ 0.99999832, -0.00183542],
                            [ 0.99417723, -0.1077573 ],
                            [ 0.99315058,  0.11684146],
                            [ 0.88772126,  0.46038133]])
        offsets = np.array([ 0.08616784, -0.62845785,  0.23081011,  0.24544835, -0.06085217, -0.82408609, -0.99721142, -0.9952619 , -0.94139908, -1.05968658, -1.18963349])
        self.assertTrue(len(E) == len(N))
        self.assertTrue(len(E) == len(O))
        self.assertTrue(E.ndim == N.ndim)
        self.assertTrue(O.ndim == 1)
        self.assertTrue(V == vert)
        self.assertTrue(E == edge)
        self.assertTrue(N == normals)
        self.assertTrue(E == offsets)
