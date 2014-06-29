from __future__ import division
import copy

import numpy as np
import scipy

from ._ds_quantize import ds_quantize
from ._evalTF import evalTF
from ._utils import diagonal_indices, _is_zpk

def simulateQDSM(u, arg2, nlev=2, x0=0):
    """v, xn, xmax, y = simulateQDSM(u, ABCD, nlev=2, x0=0)
     or
    v, xn, xmax, y = simulateQDSM(u, ntf, nlev=2, x0=0)

    Compute the output of a quadrature delta-sigma modulator with input u,
    a structure described by ABCD, an initial state x0 (default zero) and 
    a quantizer whose number of levels is specified by nlev.
    For multiple quantizers, make nlev a column vector;
    for complex quantization to a diamond lattice, multiply nlev by 1j.
     size(u) = [nu N], size(nlev) = [nq 1], size(ABCD) = [order+nq order+nq+nu]

    Alternatively, the modulator may be described by an NTF.
    The NTF is zpk object. (The STF is assumed to be 1.)
    """
    nu = u.shape[0]
    nlev = np.atleast_1d(nlev)
    nq = max(nlev.shape)
    if isinstance(arg2, scipy.signal.lti):
        k = arg2.k
        zeros = np.asarray(arg2.z)
        poles = np.asarray(arg2.p)
        form = 2
        order = max(zeros.shape)
    elif _is_zpk(arg2):
        zeros, poles, k = copy.deepcopy(arg2)
        zeros = np.asarray(zeros)
        poles = np.asarray(poles)
        form = 2
        order = max(zeros.shape)
    elif isinstance(arg2, np.ndarray): # ABCD
        if arg2.shape[1] > 2 and arg2.shape[1] == nu + arg2.shape[0]:
            # ABCD dimesions OK
            form = 1
            ABCD = arg2
            order = ABCD.shape[0] - nq
        else:
            raise ValueError('The ABCD argument does not have proper dimensions.')
    else:
        raise TypeError('The second argument is neither an ABCD matrix nor an NTF.')

    if x0 is None:
        x0 = np.zeros(shape=(order, ), dtype='float64')
    if form == 1:
        A, B, C, D = partitionABCD(ABCD, nq + nu)
        D1 = D[:, :nu]
    else:
        # Create a FF realization of 1-1/H. 
        # Note that MATLAB's zp2ss and canon functions don't work for complex TFs.
        A = np.zeros(shape=(order, order), dtype='float64')
        B2 = np.vstack((np.atleast_2d(1), np.zeros(shape=(order - 1, 1), dtype='float64')))
        diag = diagonal_indices(A, 0)
        A[diag] = zeros
        subdiag = diagonal_indices(A, -1)
        A[subdiag] = 1.
        # Compute C st C*inv(zI-A)*B = 1-1/H(z);
        w = 2*np.pi*np.random.rand(1, 2*order)
        desired = 1 - 1.0 / evalTF((zeros, poles, k), np.exp(1j*w))
        # suppress warnings about complex TFs
        sys = scipy.signal.lti(A, B2, np.eye(order), np.zeros((order, 1)))
        sysresp = reshape(freqresp(sys, w), order, max(w.shape))
        C = desired / sysresp
        # !!!! Assume stf=1
        B1 = -B2
        B = np.array([B1, B2]).reshape(1, -1)
        D1 = 1
    N = max(u.shape)
    v = np.zeros(shape=(nq, N), dtype='float64')
    y = np.zeros(shape=(nq, N), dtype='float64')
    # Need to store the state information
    xn = np.zeros(shape=(order, N), dtype='float64')
    # Need to keep track of the state maxima
    xmax = abs(x0)
    for i in range(N):
        y[:, i] = np.dot(C, x0) + np.dot(D1, u[:, i])
        v[:, i] = ds_qquantize(y[:, i], nlev)
        x0 = np.dot(A, x0) + np.dot(B, np.hstack((u[:, i], v[:, i])))
        # Save the next state
        xn[:, i] = x0
        # Keep track of the state maxima
        xmax = np.max(abs(x0), xmax)
    return v, xn, xmax, y

def ds_qquantize(y, n=2):
    """Quadrature quantization
    """
    if np.isreal(n):
        v = ds_quantize(np.real(y), n) + 1j*ds_quantize(np.imag(y), n)
    else:
        ytmp = np.array([np.real(y) + np.imag(y), np.real(y) - np.imag(y)]).reshape(1, -1)
        v = np.dot(ds_quantize(ytmp, abs(n)), np.array([1 + 1j, 1 - 1j]).reshape(1, -1)) / 2
    return v

def test_simulateQDSM():
    simulateQDSM(np.zeros(100), ((1.5,1.3), (0.2,0.2), 1))
