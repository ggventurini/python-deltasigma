# -*- coding: utf-8 -*-
# _realizeQNTF.py
# Module providing the realizeQNTF function
# Copyright 2015 Giuseppe Venturini
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

"""Module providing the realizeQNTF() function
"""

from __future__ import division, print_function

import numpy as np

from numpy.random import rand

from ._evalTF import evalTF
from ._partitionABCD import partitionABCD
from ._utils import _get_zpk, diagonal_indices

def realizeQNTF(ntf, form='FB', rot=False, bn=0.):
    """Convert a quadrature NTF into an ABCD matrix

    The basic idea is to equate the value of the loop filter at a set of
    points in the z-plane to the values of :math:`L_1 = 1-1/H` at those points.

    The order of the zeros is used when mapping the NTF onto the chosen
    topology.

    Supported structures are:

    * 'FB': Feedback
    * 'PFB': Parallel feedback
    * 'FF': Feedforward (bn is the coefficient of the aux DAC)
    * 'PFF': Parallel feedforward

    Not currently supported - feel free to send in a patch:

    * 'FBD': FB with delaying quantizer
    * 'FFD': FF with delaying quantizer

    **Parameters:**

    ntf : lti object or supported description
        The NTF to be converted.
    form : str, optional
        The form to be used. See above for the currently supported
        topologies. Defaults to ``'FB'``.
    rot : bool, optional
        Set ``rot`` to ``True`` to rotate the states to make as
        many coefficients as possible real.
    bn : float, optional
        Coefficient of the auxiliary DAC, to be specified for a
        'FF' form. Defaults to ``0.0``.

    **Returns:**

    ABCD : ndarray
        ABCD realization of the requested type.

    """
    #Code common to all forms
    ntf_z, ntf_p, _ = _get_zpk(ntf)
    order = ntf_p.shape[0]
    form = form.upper()

    A = np.diag(ntf_z)
    subdiag = diagonal_indices(A, -1)
    A[subdiag] = 1.
    if form == 'PFB' or form == 'PFF':
        # partition the NTF zeros into in-band and image-band
        for _ in range(2):
            # if the zeros are not correctly sorted, we sort them here
            fz = np.angle(ntf_z)/(2*np.pi)
            f0 = fz[0]
            inband_zeros = np.abs(fz - f0) < np.abs(fz + f0)
            ntf_z, inband_zeros = _move_zeros_to_the_end(ntf_z, inband_zeros)
        n_in = np.sum(inband_zeros)
        imageband_zeros = np.logical_not(inband_zeros)
        n_im = np.sum(imageband_zeros)
        if np.any(np.logical_not(imageband_zeros[n_in:])):
            # this should never happen.
            raise ValueError('Please put the image-band zeros at the end of' +
                             ' the ntf zeros')
        if n_im > 0:
            A[n_in, n_in-1] = 0.
    D = np.zeros(shape=(1, 2), dtype='complex128')

    # Find a set of points in the z-plane that are not close to zeros of H
    zSet = np.zeros((2*order, ), dtype='complex128')
    for i in range(2*order):
        z = 2*rand(1) - 1 + 1j*(2*rand(1) - 1)
        while np.any(np.abs(z - ntf_z) < 0.1) or np.any(np.abs(z - zSet)) < 0.1:
            z = 2*rand(1) - 1 + 1j*(2*rand(1) - 1)
        zSet[i] = z[0]
    # Evaluate L1 = 1-1/H at these points
    L1 = 1.0 - 1.0/evalTF(ntf, zSet)

    if form == 'FB':
        B = np.zeros(shape=(order, 2), dtype='complex128')
        C = np.hstack((np.zeros((1, order-1), dtype='complex128'),
                       np.atleast_2d(1)))
        # Compute F = C * inv(zI-A) at each z in zSet
        F = np.zeros((zSet.shape[0], order), dtype='complex128')
        I = np.eye(order)
        for i in range(zSet.shape[0]):
            F[i, :] = np.dot(C, np.linalg.inv(zSet[i]*I - A))
        B[:, 1] = np.linalg.lstsq(F, L1)[0]
        if rot == True:
            ABCD = np.vstack((np.hstack((A, B[:, 1].reshape((-1, 1)))),
                              np.hstack((C, np.atleast_2d(0)))))
            for i in range(order):
                phi = np.angle(ABCD[i, -1])
                ABCD[i, :] = ABCD[i, :]*np.exp(-1j*phi)
                ABCD[:, i] = ABCD[:, i]*np.exp(+1j*phi)
            A, B2, C, _ = partitionABCD(ABCD)
            B[:, 1] = np.squeeze(B2)
        B[0, 0] = np.abs(B[0, 1])
    elif form == 'PFB':
        B = np.zeros((order, 2), dtype='complex128')
        C = np.hstack((np.zeros((1, n_in-1)),
                       np.ones((1, 1)),
                       np.zeros((1, n_im-1)),
                       np.ones((1, 1))))
        # Compute F = C * inv(zI-A) at each z in zSet
        F = np.zeros((zSet.shape[0], order), dtype='complex128')
        I = np.eye(order)
        for i in range(zSet.shape[0]):
            F[i, :] = np.dot(C, np.linalg.inv(zSet[i]*I - A))
        B[:, 1] = np.linalg.lstsq(F, L1)[0]
        if rot == True:
            ABCD = np.vstack((np.hstack((A, B[:, 1].reshape((-1, 1)))),
                              np.hstack((C, np.atleast_2d(0)))))
            for i in range(order):
                phi = np.angle(ABCD[i, -1])
                ABCD[i, :] = ABCD[i, :]*np.exp(-1j*phi)
                ABCD[:, i] = ABCD[:, i]*np.exp(+1j*phi)
            A, B2, C, _ = partitionABCD(ABCD)
            B[:, 1] = np.squeeze(B2)
        B[0, 0] = np.abs(B[0, 1])
    elif form == 'FF':
        B0 = np.vstack((np.atleast_2d(1),
                        np.zeros((order-1, 1), dtype='complex128')))
        B = B0.copy()
        B[-1, 0] = bn
        # Compute F = inv(zI-A)*B at each z in zSet
        F = np.zeros((order, zSet.shape[0]), dtype='complex128')
        I = np.eye(order)
        for i in range(zSet.shape[0]):
            F[:, i] = np.squeeze(np.dot(np.linalg.inv(zSet[i]*I - A), B))
        L1 = L1.reshape((-1, 1))
        C = np.linalg.lstsq(F.T, L1)[0].T
        if rot == True:
            ABCD = np.vstack((np.hstack((A, B)),
                              np.hstack((C, np.atleast_2d(0)))))
            for i in range(1, order-1):
                phi = np.angle(ABCD[-1, i])
                ABCD[i, :] = ABCD[i, :]*np.exp(+1j*phi)
                ABCD[:, i] = ABCD[:, i]*np.exp(-1j*phi)
            A, B, C, _ = partitionABCD(ABCD)
        B = np.hstack((-B0, B))
    elif form == 'PFF':
        B0 = np.vstack((np.atleast_2d(1),
                        np.zeros((order-1, 1), dtype='complex128')))
        B = B0.copy()
        B[n_in] = 1.
        # Compute F = inv(zI-A)*B at each z in zSet
        F = np.zeros((order, zSet.shape[0]), dtype='complex128')
        I = np.eye(order)
        for i in range(zSet.shape[0]):
            F[:, i] = np.squeeze(np.dot(np.linalg.inv(zSet[i]*I - A), B))
        L1 = L1.reshape((-1, 1))
        C = np.linalg.lstsq(F.T, L1)[0].T
        C = C[:, ::-1]
        if rot == True:
            ABCD = np.vstack((np.hstack((A, B)),
                              np.hstack((C, np.zeros((1, 1))))))
            for i in range(1, order-1):
                phi = np.angle(ABCD[-1, i])
                ABCD[i, :] = ABCD[i, :]*np.exp(+1j*phi)
                ABCD[:, i] = ABCD[:, i]*np.exp(-1j*phi)
            A, B, C, _ = partitionABCD(ABCD)
        B = np.hstack((B0, B))
    else:
        raise ValueError('Sorry, form %s is not supported.' % form)
    ABCD = np.vstack((np.hstack((A, B)), np.hstack((C, D))))
    return ABCD

def _move_zeros_to_the_end(ntf_z, inband_zeros):
    """RealizeQNTF requires the imageband zeros to be put at the end
    of the list of zeros.

    This function ensures exactly that.
    """
    to_the_end = [i for i, iz in enumerate(inband_zeros) if not iz]
    for i, orig in enumerate(to_the_end):
        ntf_z = np.append(ntf_z, ntf_z[orig - i])
        ntf_z = np.delete(ntf_z, orig - i)
        inband_zeros = np.append(inband_zeros, inband_zeros[orig - i])
        inband_zeros = np.delete(inband_zeros, orig - i)
    return ntf_z, inband_zeros
