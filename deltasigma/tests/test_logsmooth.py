# -*- coding: utf-8 -*-
# test_logsmooth.py
# This module provides the tests for the logsmooth function.
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

"""This module provides the test class for the logsmooth() function.
"""

from __future__ import division

import unittest
import numpy as np
import pkg_resources, scipy.io

from numpy.fft import fft

from deltasigma import ds_f1f2
from deltasigma import ds_hann
from deltasigma import logsmooth
from deltasigma import simulateDSM
from deltasigma import synthesizeNTF
from deltasigma import undbv

class TestLogsmooth(unittest.TestCase):
    """Test class for logsmooth()"""

    def setUp(self):
        pass

    def test_logsmooth(self):
        """Test function for logsmooth()"""
        f0 = 1./8
        OSR = 64
        order = 8
        N = 8192
        H = synthesizeNTF(order, OSR, 1, 1.5, f0)
        # fB = int(np.ceil(N/(2. * OSR)))
        quadrature = False
        M = 1
        f1, f2 = ds_f1f2(OSR, f0, quadrature)
        Amp = undbv(-3)
        f = .3
        fin = np.round(((1 - f)/2*f1 + (f + 1)/2 * f2) * N)
        t = np.arange(0, N).reshape((1, -1))
        u = Amp * M * np.cos((2*np.pi/N)*fin*t)
        v, _, _, _ = simulateDSM(u, H, M + 1)

        window = ds_hann(N)
        # NBW = 1.5/N
        spec0 = fft(v * window) / (M*N/4)

        ### -1 is really important THIS IS NOT MATLAB!!
        fl, pl = logsmooth(spec0, fin - 1)

        fname = pkg_resources.resource_filename(__name__,
                                                "test_data/test_logsmooth.mat")
        data = scipy.io.loadmat(fname)
        self.assertTrue(np.allclose(fl, data['fl'], atol=1e-8, rtol=1e-5))
        self.assertTrue(np.allclose(pl, data['pl'], atol=1e-8, rtol=1e-5))


