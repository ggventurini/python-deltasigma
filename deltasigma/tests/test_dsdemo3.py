from __future__ import division
import unittest
import numpy as np
import deltasigma as ds
import pkg_resources
import scipy as sp
import scipy.io


## Modulator realization and dynamic range scaling - # demo3
# In this ipython notebook, the following is demonstrated:
# 
#  * A 5th order delta sigma modulator is synthesized, with optimized zeros and an OSR equal to 42.
# 
#  * We then convert the synthesized NTF into `a`, `g`, `b`, `c` coefficients for the `CRFB` modulator structure.
#  
#  * The maxima for each state are evaluated.
#  
#  * The `ABCD` matrix is scaled so that the state maxima are less than the specified limit.
# 
#  * The state maxima are re-evaluated and limit compliance is checked.
# 
# **NOTE:** This is an ipython port of `dsdemo3.m`, from the
# **[MATLAB Delta Sigma Toolbox](http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox)**, 
# written by Richard Schreier.

class TestAxisLabels(unittest.TestCase):
    """dsdemo3: Delta sigma modulator synthesis"""
    def setUp(self):
        fname = pkg_resources.resource_filename(__name__, "test_data/test_dsdemo3.mat")
        self.data = sp.io.loadmat(fname)
    def test_dsdemo3(self):
        """dsdemo3 test: Delta sigma modulator synthesis"""
        order = 5
        R = 42
        opt = 1
        H = ds.synthesizeNTF(order, R, opt)

        # ## Evaluation of the coefficients for a CRFB topology
        a, g, b, c = ds.realizeNTF(H)

        # Use a single feed-in for the input
        # Lets check that stuffABCD understands that if b is scalar it means 1 feed-in.
        ABCD = ds.stuffABCD(a, g, b[0], c)
        # for passing the assertion below, we need b to have the trailing zeros
        b = np.concatenate((np.atleast_1d(b[0]), 
                            np.zeros((b.shape[0] - 1,))))
        u = np.linspace(0, 0.6, 30)
        N = 1e4
        T = np.ones((1, N))
        maxima = np.zeros((order, len(u)))
        for i in range(len(u)):
            ui = u[i]
            v, xn, xmax, _ = ds.simulateDSM(ui*T, ABCD);
            maxima[:, i] = np.squeeze(xmax)
            if any(xmax > 1e2):
                umax = ui
                u = u[:i+1]
                maxima = maxima[:, :i]
                break
        # save the maxima
        prescale_maxima = np.copy(maxima)

        # ## Scaled modulator
        # ### Calculate the scaled coefficients
        ABCDs, umax, _ = ds.scaleABCD(ABCD, N_sim=1e5)
        as_, gs, bs, cs = ds.mapABCD(ABCDs)
        # ### Calculate the state maxima
        u = np.linspace(0, umax, 30)
        N = 1e4
        T = np.ones((N,))
        maxima = np.zeros((order, len(u)))
        for i in range(len(u)):
            ui = u[i]
            v, xn, xmax, _ = ds.simulateDSM(ui*T, ABCDs)
            maxima[:, i] = xmax.squeeze()
            if any(xmax > 1e2):
                umax = ui;
                u = u[:i]
                maxima = maxima[:, :i]
                break

        self.assertTrue(np.allclose(ABCD, self.data['ABCD']))
        self.assertTrue(np.allclose(ABCDs, self.data['ABCDs'], atol=1e-2, rtol=5e-2))
        self.assertTrue(np.allclose(a, self.data['a'], atol=1e-5, rtol=1e-3))
        self.assertTrue(np.allclose(b, self.data['b'], atol=1e-5, rtol=1e-3))
        self.assertTrue(np.allclose(g, self.data['g'], atol=1e-5, rtol=1e-3))
        self.assertTrue(np.allclose(c, self.data['c'], atol=1e-5, rtol=1e-3))
        self.assertTrue(np.allclose(as_, self.data['as'], atol=1e-2, rtol=5e-3))
        self.assertTrue(np.allclose(bs, self.data['bs'], atol=5e-3, rtol=5e-3))
        self.assertTrue(np.allclose(gs, self.data['gs'], atol=5e-3, rtol=5e-3))
        self.assertTrue(np.allclose(cs, self.data['cs'], atol=3e-2, rtol=3e-2))
        self.assertTrue(np.allclose(umax, self.data['umax'], atol=5e-3, rtol=5e-3))
        self.assertTrue(np.allclose(maxima, self.data['maxima'], atol=2e-2, rtol=5e-2))

