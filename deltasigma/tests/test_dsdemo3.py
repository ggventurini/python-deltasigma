from __future__ import division
from ._synthesizeNTF import synthesizeNTF
from ._realizeNTF import realizeNTF
from ._stuffABCD import stuffABCD
from ._simulateDSM import simulateDSM
from ._scaleABCD import scaleABCD
from ._mapABCD import mapABCD
import pkg_resources
import scipy.io
import numpy as np

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
# **NOTE:** This is an ipython port of `dsdemo3.m`, from the **[MATLAB Delta Sigma Toolbox](http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox)**, written by Richard Schreier.

def test_dsdemo3():
    # ## Delta sigma modulator synthesis
    order = 5
    R = 42
    opt = 1
    H = synthesizeNTF(order, R, opt)

    # ## Evaluation of the coefficients for a CRFB topology
    a, g, b, c = realizeNTF(H)

    # Use a single feed-in for the input
    # Lets check that stuffABCD understands that if b is scalar it means 1 feed-in.
    ABCD = stuffABCD(a, g, b[0], c)
    # for passing the assertion below, we need b to have the trailing zeros
    b = np.concatenate((np.atleast_1d(b[0]), 
                        np.zeros((b.shape[0] - 1,))))
    u = np.linspace(0, 0.6, 30)
    N = 1e4
    T = np.ones((1, N))
    maxima = np.zeros((order, len(u)))
    for i in range(len(u)):
        ui = u[i]
        v, xn, xmax, _ = simulateDSM(ui*T, ABCD);
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
    ABCDs, umax, _ = scaleABCD(ABCD, N_sim=1e5)
    as_, gs, bs, cs = mapABCD(ABCDs)
    # ### Calculate the state maxima
    u = np.linspace(0, umax, 30)
    N = 1e4
    T = np.ones((N,))
    maxima = np.zeros((order, len(u)))
    for i in range(len(u)):
        ui = u[i]
        v, xn, xmax, _ = simulateDSM(ui*T, ABCDs)
        maxima[:, i] = xmax.squeeze()
        if any(xmax > 1e2):
            umax = ui;
            u = u[:i]
            maxima = maxima[:, :i]
            break

    fname = pkg_resources.resource_filename(__name__, "test_data/test_dsdemo3.mat")
    data = scipy.io.loadmat(fname)

    assert np.allclose(ABCD, data['ABCD'])
    assert np.allclose(ABCDs, data['ABCDs'], atol=1e-2, rtol=5e-2)
    assert np.allclose(a, data['a'], atol=1e-5, rtol=1e-3)
    assert np.allclose(b, data['b'], atol=1e-5, rtol=1e-3)
    assert np.allclose(g, data['g'], atol=1e-5, rtol=1e-3)
    assert np.allclose(c, data['c'], atol=1e-5, rtol=1e-3)
    assert np.allclose(as_, data['as'], atol=1e-2, rtol=5e-3)
    assert np.allclose(bs, data['bs'], atol=5e-3, rtol=5e-3)
    assert np.allclose(gs, data['gs'], atol=5e-3, rtol=5e-3)
    assert np.allclose(cs, data['cs'], atol=3e-2, rtol=3e-2)
    assert np.allclose(umax, data['umax'], atol=5e-3, rtol=5e-3)
    assert np.allclose(maxima, data['maxima'], atol=2e-2, rtol=5e-2)
