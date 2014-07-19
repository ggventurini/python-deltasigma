from __future__ import division
import unittest
import numpy as np
import deltasigma as ds
import scipy.io
import pkg_resources
from deltasigma._utils import mround


class TestBplogsmooth(unittest.TestCase):
    """Test class for bplogsmooth()"""

    def setUp(self):
        fname = pkg_resources.resource_filename(
            __name__, "test_data/test_bplogsmooth.mat")
        self.data = scipy.io.loadmat(fname)
        f0 = 1./8
        OSR = 64
        order = 8
        N = 8192
        H = ds.synthesizeNTF(order, OSR, 1, 1.5, f0)
        fB = int(np.ceil(N/(2. * OSR)))
        ftest = int(mround(f0*N + 1./3*fB))
        u = 0.5*np.sin(2*np.pi*ftest/N*np.arange(N))
        v, xn, xmax, y = ds.simulateDSM(u, H)
        spec = np.fft.fft(v*ds.ds_hann(N))/(N/4)
        X = spec[:N/2 + 1]
        self.f, self.p = ds.bplogsmooth(X, ftest, f0)

    def test_one(self):
        """Test function 1/2 for bplogsmooth()"""
        data = self.data
        f = self.f
        self.assertTrue(np.allclose(f, data['f'], atol=1e-9, rtol=1e-5))

    def test_two(self):
        """Test function 2/2 for bplogsmooth()"""
        data = self.data
        p = self.p
        self.assertTrue(np.allclose(p, data['p'], atol=1e-9, rtol=1e-5))
