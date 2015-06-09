# -*- coding: utf-8 -*-
# _calculateTF.py
# Module providing the calculateTF function
# Copyright 2013 Giuseppe Venturini
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

"""Module providing the calculateTF() function
"""

from warnings import warn

import numpy as np
from scipy.signal import lti, ss2zpk

from ._constants import eps
from ._partitionABCD import partitionABCD
from ._utils import carray, minreal


def calculateTF(ABCD, k=1.):
    """Calculate the NTF and STF of a delta-sigma modulator.

    The calculation is performed for a given loop filter
    ABCD matrix, assuming a quantizer gain of ``k``.

    **Parameters:**

    ABCD : array_like,
        The ABCD matrix that describes the system.

    k : float or ndarray-like, optional
        The quantizer gains. If only one quantizer is present, it may be set
        to a float, corresponding to the quantizer gain. If multiple quantizers
        are present, a list should be used, with quantizer gains ordered
        according to the order in which the quantizer inputs apperar in the
        ``C`` and ``D`` submatrices. If not specified, a default of one quatizer
        with gain ``1.`` is assumed.

    **Returns:**

    (NTF, STF) : a tuple of two LTI objects (or of two lists of LTI objects).

    If the system has multiple quantizers, multiple STFs and NTFs will be
    returned.

    In that case:

    * ``STF[i]`` is the STF from ``u`` to output number ``i``.
    * ``NTF[i, j]`` is the NTF from the quantization noise of the quantizer
      number ``j`` to output number ``i``.

    **Note:**

    Setting ``k`` to a list is unsupported in the MATLAB code (last checked
    Nov. 2014).

    **Example:**

    Realize a fifth-order modulator with the cascade-of-resonators structure,
    feedback form. Calculate the ABCD matrix of the loop filter and verify
    that the NTF and STF are correct.

    .. code-block:: python

        from deltasigma import *
        H = synthesizeNTF(5, 32, 1)
        a, g, b, c = realizeNTF(H)
        ABCD = stuffABCD(a,g,b,c)
        ntf, stf = calculateTF(ABCD)

    From which we get:

    ``H``::

             (z -1) (z^2 -1.997z +1) (z^2 -1.992z +0.9999)
        --------------------------------------------------------
         (z -0.7778) (z^2 -1.796z +0.8549) (z^2 -1.613z +0.665)

    coefficients::

        a: 0.0007, 0.0084, 0.055, 0.2443, 0.5579
        g: 0.0028, 0.0079
        b: 0.0007, 0.0084, 0.055, 0.2443, 0.5579, 1.0
        c: 1.0, 1.0, 1.0, 1.0, 1.0

    ABCD matrix::

        [[  1.00000000e+00   0.00000000e+00   0.00000000e+00   0.00000000e+00
            0.00000000e+00   6.75559806e-04  -6.75559806e-04]
         [  1.00000000e+00   1.00000000e+00  -2.79396240e-03   0.00000000e+00
            0.00000000e+00   8.37752565e-03  -8.37752565e-03]
         [  1.00000000e+00   1.00000000e+00   9.97206038e-01   0.00000000e+00
            0.00000000e+00   6.33294166e-02  -6.33294166e-02]
         [  0.00000000e+00   0.00000000e+00   1.00000000e+00   1.00000000e+00
           -7.90937431e-03   2.44344030e-01  -2.44344030e-01]
         [  0.00000000e+00   0.00000000e+00   1.00000000e+00   1.00000000e+00
            9.92090626e-01   8.02273699e-01  -8.02273699e-01]
         [  0.00000000e+00   0.00000000e+00   0.00000000e+00   0.00000000e+00
            1.00000000e+00   1.00000000e+00   0.00000000e+00]]

    NTF::

             (z -1) (z^2 -1.997z +1) (z^2 -1.992z +0.9999)
        --------------------------------------------------------
         (z -0.7778) (z^2 -1.796z +0.8549) (z^2 -1.613z +0.665)

    STF::

        1

    """

    nq = len(k) if type(k) in (tuple, list) else 1
    A, B, C, D = partitionABCD(ABCD, m=nq+1, r=nq)
    k = carray(k)
    diagk = np.atleast_2d(np.diag(k))

    B1 = B[:, 0]
    B2 = B[:, 1:]
    B1 = B1.reshape((B1.shape[0], 1)) if len(B1.shape) == 1 else B1
    B2 = B2.reshape((B2.shape[0], 1)) if len(B2.shape) == 1 else B2

    # In a single-quantizer system, D2 should be all zeros, as any
    # non-zero term in D2 would imply we got a delay-free loop.
    # The original MATLAB code implies so.
    # The trouble arises when we adapt and extend the code to calculate
    # the transfer functions for multiple-quantizer systems. There the
    # off-diagonal terms of D2 may very well be non-zero, therefore
    # in the following we consider D2.
    # this means that if you supply an ABCD matrix, with one quantizer
    # and an (erroneously) zero-delay loop, the MATLAB toolbox will
    # disregard D2 and pretend it's zero and the loop is correctly
    # delayed, spitting out the corresponding TF.
    # Instead here we print out a warning and process the ABCD matrix
    # with the user-supplied, non-zero D2 matrix.
    # The resulting TFs will obviously be different.

    D1 = D[:, 0]
    D2 = D[:, 1:]
    D1 = D1.reshape((D1.shape[0], 1)) if len(D1.shape) == 1 else D1
    D2 = D2.reshape((D2.shape[0], 1)) if len(D2.shape) == 1 else D2

    # WARN DELAY FREE LOOPS
    if np.diag(D2).any():
        warn("Delay free loop detected! D2 diag: %s", str(np.diag(D2)))

    # Find the noise transfer function by forming the closed-loop
    # system (sys_cl) in state-space form.
    Ct = np.linalg.inv(np.eye(nq) - D2*diagk)
    Acl = A + np.dot(B2, np.dot(Ct, np.dot(diagk, C)))
    Bcl = np.hstack((B1 + np.dot(B2, np.dot(Ct, np.dot(diagk, D1))),
                     np.dot(B2, Ct)))
    Ccl = np.dot(Ct, np.dot(diagk, C))
    Dcl = np.dot(Ct, np.hstack((np.dot(diagk, D1), np.eye(nq))))
    tol = min(1e-3, max(1e-6, eps**(1/ABCD.shape[0])))
    ntfs = np.empty((nq, Dcl.shape[0]), dtype=np.object_)
    stfs = np.empty((Dcl.shape[0],), dtype=np.object_)

    # sweep the outputs 'cause scipy is silly but we love it anyhow.
    for i in range(Dcl.shape[0]):
        # input #0 is the signal
        # inputs #1,... are quantization noise
        stf_z, stf_p, stf_k  = ss2zpk(Acl, Bcl, Ccl[i, :], Dcl[i, :], input=0)
        stf = lti(stf_z, stf_p, stf_k)
        for j in range(nq):
            ntf_z, ntf_p, ntf_k = ss2zpk(Acl, Bcl, Ccl[i, :], Dcl[i, :], input=j+1)
            ntf = lti(ntf_z, ntf_p, ntf_k)
            stf_min, ntf_min = minreal((stf, ntf), tol)
            ntfs[i, j] = ntf_min
        stfs[i] = stf_min

    # if we have one stf and one ntf, then just return those in a list
    if ntfs.shape == (1, 1):
        return [ntfs[0, 0], stfs[0]]
    return ntfs, stfs

