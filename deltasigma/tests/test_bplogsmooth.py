import unittest
import numpy as np
import deltasigma as ds
import scipy.io
import pkg_resources
from deltasigma._ds_hann import ds_hann
from deltasigma._simulateDSM import simulateDSM
from deltasigma._synthesizeNTF import synthesizeNTF
from deltasigma._utils import mround


class TestBplogsmooth(unittest.TestCase):
    """Test function for bplogsmooth()"""

    def setUp(self):
        fname = pkg_resources.resource_filename(
            __name__, "test_data/test_bplogsmooth.mat")
        self.data = scipy.io.loadmat(fname)
        f0 = 1./8
        OSR = 64
        order = 8
        N = 8192
        H = synthesizeNTF(order, OSR, 1, 1.5, f0)
        fB = int(np.ceil(N/(2. * OSR)))
        ftest = int(mround(f0*N + 1./3*fB))
        u = 0.5*np.sin(2*np.pi*ftest/N*np.arange(N))
        v, xn, xmax, y = simulateDSM(u, H)
        spec = np.fft.fft(v*ds_hann(N))/(N/4)
        X = spec[:N/2 + 1]
        self.f, self.p = ds.bplogsmooth(X, ftest, f0)

    def test_one(self):
        data = self.data
        f = self.f
        self.assertTrue(np.allclose(f, data['f'], atol=1e-9, rtol=1e-5))

    def test_two(self):
        data = self.data
        p = self.p
        self.assertTrue(np.allclose(p, data['p'], atol=1e-9, rtol=1e-5))
