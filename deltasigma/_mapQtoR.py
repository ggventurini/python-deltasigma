# -*- coding: utf-8 -*-
# _mapQtoR.py
# Module providing the mapQtoR function
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

"""Module providing the mapQtoR() function
"""

from __future__ import division
import numpy as np

def mapQtoR(ABCD):
    """Map a quadrature ABCD matrix to a real one.

    Each element z in ABCD is replaced by a 2x2 matrix in ``ABCDr``, the
    return value.

    Specifically:

    .. math::

        z \\rightarrow
         \\begin{bmatrix}
          x & -y \\\\
          y & x \\\\
         \\end{bmatrix}
        \\mathrm{where}\\ x = Re(z)\\ \\mathrm{and}\\ y = Im(z)

    The non-quadrature topology can be simulated with :func:`simulateDSM`.::

        import numpy as np
        from deltasigma import *
        nlev = 9
        f0 = 1./16.
        osr = 32
        M = nlev - 1
        ntf = synthesizeQNTF(4, osr, f0, -50, -10)
        N = 64*osr
        f = int(np.round((f0 + 0.3*0.5/osr)*N)/N)
        u = 0.5*M*np.exp(2j*np.pi*f*np.arange(N))
        # Instead of calling simulateQDSM
        # v = simulateQDSM(u, ntf, nlev)
        # it's faster to run:
        ABCD = realizeQNTF(ntf, 'FF')
        ABCDr = mapQtoR(ABCD)
        ur = np.vstack((np.real(u), np.imag(u)))
        vr = simulateDSM(ur, ABCDr, nlev*np.array([[1],[1]]))
        v = vr[0,:] + 1j*vr[1, :]

    Notice the example above requires the function :func:`synthesizeQNTF`,
    which is not part of the current release of python-deltasigma.

    """
    A = np.zeros((2*ABCD.shape[0], 2*ABCD.shape[1]))
    A[::2, ::2] = np.real(ABCD)
    A[1::2, 1::2] = np.real(ABCD)
    A[::2, 1::2] = -np.imag(ABCD)
    A[1::2, ::2] = +np.imag(ABCD)
    return A

