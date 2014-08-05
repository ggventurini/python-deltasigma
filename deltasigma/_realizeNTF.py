# -*- coding: utf-8 -*-
# _realizeNTF.py
# Module providing realizeNTF
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

"""Module providing realizeNTF()
"""

from __future__ import division, print_function
from warnings import warn

import numpy as np

from ._calculateTF import calculateTF
from ._evalRPoly import evalRPoly
from ._evalTF import evalTF
from ._stuffABCD import stuffABCD
from ._utils import carray, cplxpair, _get_zpk


def realizeNTF(ntf, form='CRFB', stf=None):
    """Convert an NTF into coefficients for the desired structure.

    **Parameters:**

    ntf : a zpk transfer function
        The modulator noise transfer function (NTF)

    form : string
        A structure identifier.

    Supported structures are:

    * *"CRFB"*: Cascade of resonators, feedback form.

    * *"CRFF"*: Cascade of resonators, feedforward form.

    * *"CIFB"*: Cascade of integrators, feedback form.

    * *"CIFF"*: Cascade of integrators, feedforward form.

    * *"CRFBD"*: CRFB with delaying quantizer.

    * *"CRFFD"*: CRFF with delaying quantizer.

    * *"PFF"*: Parallel feed-forward.

    * *"Stratos"*: A CIFF-like structure with non-delaying resonator feedbacks,
      contributed to the MATLAB Delta Sigma toolbox in 2007 by Jeff Gealow.

    See the accompanying documentation (:ref:`topologies-diagrams`) for
    block diagrams of each structure.

    .. note:

        The order of the NTF zeros must be (real, complex conj. pairs).
        The order of the zeros is used when mapping the NTF onto the chosen
        topology.

    stf : a zpk transfer function, optional
        the Signal Transfer Function

    **Returns:**

    a, g, b, c : tuple of ndarrays
        the coefficients for the desired structure

    **Example:**

    Determine the coefficients for a 5th-order modulator with the
    cascade-of-resonators structure, feedback (CRFB) form.::

        from deltasigma import synthesizeNTF, realizeNTF
        H = synthesizeNTF(5, 32, 1)
        a, g, b, c = realizeNTF(H,'CRFB')

    Returns the values::

        a: 0.0007, 0.0084, 0.055, 0.2443, 0.5579
        g: 0.0028, 0.0079
        b: 0.0007, 0.0084, 0.055, 0.2443, 0.5579, 1.0
        c: 1.0, 1.0, 1.0, 1.0, 1.0

    """

    # The basic idea is to equate the loop filter at a set of
    # points in the z-plane to L1  =  1-1/ntf at those points.

    # Code common to all functions\
    if (hasattr(ntf, 'inputs') and not ntf.inputs == 1) or \
       (hasattr(ntf, 'outputs') and not ntf.outputs == 1):
        raise TypeError("Only SISO transfer functions can be evaluated.")

    ntf_z, ntf_p, _ = _get_zpk(ntf)
    ntf_z = carray(ntf_z)
    ntf_p = carray(ntf_p)
    order = max(ntf_p.shape)
    order2 = int(np.floor(order / 2))
    odd = (order - 2 * order2 == 1)

    a = np.zeros((1, order))
    g = np.zeros((1, order2))
    b = np.zeros((1, order + 1))
    c = np.ones((1, order))
    T = np.zeros((order, order * 2), dtype='complex')

    # instead of asuming enforce roots order
    ntf_z = cplxpair(ntf_z)[::-1]

    N = 200
    min_distance = 0.09
    zSet = np.zeros((N, ), dtype='complex')
    j = 0
    for i in range(1, N + 1):
        z = 1.1 * np.exp(2j * np.pi * i / N)
        if np.all(np.abs(ntf_z - z) > min_distance):
            zSet[j] = z
            j = j + 1
    zSet = zSet[:j]
    if form == 'CRFB':
        # Find g
        # Assume the roots are ordered, real first, then cx conj. pairs
        for i in range(0, order2):
            g[0, i] = 2 * (1 - np.real(ntf_z[2 * i + 1 * odd]))
        L1 = np.zeros((1, order * 2), dtype='complex')
        # Form the linear matrix equation a*T* = L1
        for i in range(0, 2 * order):
            z = zSet[i]
            L1[0, i] = 1. - evalRPoly(ntf_p, z) / evalRPoly(ntf_z, z)
            Dfactor = (z - 1) / z
            product = 1
            for j in range(order, 0 + 1 * odd, -2):
                product = z / evalRPoly(ntf_z[j - 2:j], z) * product
                T[j - 1, i] = product * Dfactor
                T[j - 2, i] = product
            if odd:
                T[0, i] = product / (z - 1)
        a = -np.real(np.linalg.lstsq(T.T, L1.T)[0]).T
        if stf is None:
            b[0,:order] = a[0, :]
            b[0, order] = 1
    elif form == 'CRFF':
        # Find g
        # Assume the roots are ordered, real first, then cx conj. pairs
        for i in range(0, order2):
            g[0, i] = 2 * (1 - np.real(ntf_z[2 * i + 1 * odd]))
        L1 = np.zeros((1, order * 2), dtype='complex')
        # Form the linear matrix equation a*T*=L1
        for i in range(0, 2 * order):
            z = zSet[i]
            # L1(z) = 1-1/H(z)
            L1[0, i] = 1. - evalRPoly(ntf_p, z) / evalRPoly(ntf_z, z)
            if odd:
                Dfactor = z - 1
                product = 1 / Dfactor
                T[0, i] = product
            else:
                Dfactor = (z - 1) / z
                product = 1.
            for j in range(0 + 1 * odd, order, 2):
                product = z / evalRPoly(ntf_z[j:j + 2], z) * product
                T[j, i] = product * Dfactor
                T[j + 1, i] = product
        a = -np.real(np.linalg.lstsq(T.T, L1.T)[0]).T
        if stf is None:
            b = np.hstack((np.atleast_2d(1),
                           np.zeros((1, order - 1)),
                           np.atleast_2d(1)
                         ))
    elif form == 'CIFB':
        if any(abs(np.real(ntf_z) - 1) > 0.001):
            warn("The ntf's zeros have had their real parts set to one.")
        ntf_z = 1 + 1j * np.imag(ntf_z)
        for i in range(order2):
            g[0, i] = np.imag(ntf_z[2 * i + 1 * odd]) ** 2
        L1 = np.zeros((1, order * 2), dtype='complex')
        for i in range(order * 2):
            z = zSet[i]
            L1[0, i] = 1. - evalRPoly(ntf_p, z) / evalRPoly(ntf_z, z)
            Dfactor = z - 1.
            product = 1.
            for j in range(order - 1, 1 * odd, -2):
                product = product / evalRPoly(ntf_z[j - 1:j + 1], z)
                T[j, i] = product * Dfactor
                T[j - 1, i] = product
            if odd:
                T[0, i] = product / (z - 1.)
        a = -np.real(np.linalg.lstsq(T.T, L1.T)[0]).T
        if stf is None:
            b[0,:order] = a
            b[0, order] = 1.
    elif form == 'CIFF':
        if (abs(np.real(ntf_z) - 1) > 0.001).any():
            warn("The ntf's zeros have had their real parts set to one.")
        ntf_z = 1. + 1j * np.imag(ntf_z)
        for i in range(order2):
            g[0, i] = np.imag(ntf_z[2 * i + 1 * odd]) ** 2
        L1 = np.zeros((1, order * 2), dtype='complex')
        for i in range(order * 2):
            z = zSet[i]
            L1[0, i] = 1. - evalRPoly(ntf_p, z) / evalRPoly(ntf_z, z)
            Dfactor = (z - 1)
            if odd:
                product = 1. / (z - 1.)
                T[0, i] = product
            else:
                product = 1
            for j in range(0 + 1 * odd, order - 1, 2):
                product = product / evalRPoly(ntf_z[j:j + 2], z)
                T[j, i] = product * Dfactor
                T[j + 1, i] = product
        a = -np.real(np.linalg.lstsq(T.T, L1.T)[0]).T
        if stf is None:
            b = np.zeros((1, order + 1))
            b[0, 0] = 1.
            b[0, -1] = 1.
    elif form == 'CRFBD':
        for i in range(order2):
            g[0, i] = 2 * (1. - np.real(ntf_z[2 * i + 1 * odd]))
        L1 = np.zeros((1, order * 2), dtype='complex')
        for i in range(order * 2):
            z = zSet[i]
            L1[0, i] = 1. - evalRPoly(ntf_p, z) / evalRPoly(ntf_z, z)
            Dfactor = (z - 1.)
            product = 1. / z
            for j in range(order - 1, 0 + 1 * odd, -2):
                product = z / evalRPoly(ntf_z[j - 1:j + 1], z) * product
                T[j, i] = product * Dfactor
                T[j - 1, i] = product
            if odd:
                T[0, i] = product * z / (z - 1.)
        a = -np.real(np.linalg.lstsq(T.T, L1.T)[0]).T
        if stf is None:
            b = np.hstack((a, np.zeros((1, 1))))
            b[0, order] = 1.
    elif form == 'CRFFD':
        for i in range(order2):
            g[0, i] = 2 * (1. - np.real(ntf_z[2 * i + 1 * odd]))
        zL1 = np.zeros((1, order*2), dtype='complex')
        # Form the linear matrix equation a*T*=zL1
        for i in range(order * 2):
            z = zSet[i]
            # zL1 = z*(1-1/H(z))
            zL1[0, i] = z*(1. - 1.0/evalTF(ntf, z))
            if odd:
                Dfactor = (z - 1.) / z
                product = 1. / Dfactor
                T[0, i] = product
            else:
                Dfactor = z - 1
                product = 1
            for j in range(0 + 1 * odd, order, 2):
                product = z / evalRPoly(ntf_z[j:j + 2], z) * product
                T[j, i] = product * Dfactor
                T[j + 1, i] = product
        a = -np.real(np.linalg.lstsq(T.T, zL1.T)[0]).T
        if stf is None:
            b = np.hstack((np.atleast_2d(1),
                           np.zeros((1, order - 1)),
                           np.atleast_2d(1))
                         )
    elif form == 'PFF':
        warn("Untested code accessed: realizeNTF('PFF')")
        # Find g
        # Assume the roots are ordered, real first, then cx conj. pairs
        # with the secondary zeros after the primary zeros
        for i in range(order2):
            g[0, i] = 2 * (1. - np.real(ntf_z[2 * i + 1 * odd]))
        # Find the dividing line between the zeros
        theta0 = np.abs(np.angle(ntf_z[0]))
        # !! 0.5 radians is an arbitrary separation !!
        i = np.flatnonzero(np.abs(np.abs(np.angle(ntf_z)) - theta0) > 0.5)
        if i.shape[0]:
            order_1 = i[0] - 1
        else:
            order_1 = 0
        order_2 = order - order_1
        if i.shape[0] != order_2:
            raise ValueError('For the PFF form, the NTF zeros must ' +
                             'be sorted into primary and secondary zeros')
        odd_1 = order_1 % 2
        odd_2 = order_2 % 2
        L1 = np.zeros((1, order * 2), dtype='complex')
        # Form the linear matrix equation a*T*=L1
        for i in range(order * 2):
            z = zSet[i]
            # L1(z) = 1-1/H(z)
            L1[0, i] = 1. - evalRPoly(ntf_p, z) / evalRPoly(ntf_z, z)
            if odd_1:
                Dfactor = z - 1
                product = 1. / Dfactor
                T[0, i] = product
            else:
                Dfactor = (z - 1) / z
                product = 1
            for j in range(odd_1, order_1, 2):
                product = z / evalRPoly(ntf_z[j:j + 2], z) * product
                T[j, i] = product * Dfactor
                T[j + 1, i] = product
            if odd_2:
                Dfactor = z - 1
                product = 1. / Dfactor
                T[order_1, i] = product
            else:
                Dfactor = (z - 1) / z
                product = 1.
            for j in range(order_1 + odd_2, order, 2):
                product = z / evalRPoly(ntf_z[j:j + 2], z) * product
                T[j, i] = product * Dfactor
                T[j + 1, i] = product
        a = -np.real(np.linalg.lstsq(T.T, zL1.T)[0]).T
        if stf is None:
            b = np.hstack((np.atleast_2d(1),
                           np.zeros((1, order_1 - 1)),
                           np.atleast_2d(1),
                           np.zeros(1, order_2 - 1),
                           np.atleast_2d(1)
                         ))
    elif form == 'Stratos':
        # code copied from case 'CRFF':
        # Find g
        # Assume the roots are ordered, real first, then cx conj. pairs
        for i in range(order2):
            g[0, i] = 2 * (1 - np.real(ntf_z[2 * i + 1 * odd]))
        # code copied from case 'CIFF':
        L1 = np.zeros((1, order * 2), dtype='complex')
        # Form the linear matrix equation a*T*=L1
        for i in range(order * 2):
            z = zSet[i]
            # L1(z) = 1-1/H(z)
            L1[0, i] = 1. - evalRPoly(ntf_p, z) / evalRPoly(ntf_z, z)
            Dfactor = (z - 1.)
            if odd:
                product = 1. / (z - 1.)
                T[0, i] = product
            else:
                product = 1.
            for j in range(0 + 1 * odd, order - 1, 2):
                product = product / evalRPoly(ntf_z[j:j + 2], z)
                T[j, i] = product * Dfactor
                T[j + 1, i] = product
        a = -np.real(np.linalg.lstsq(T.T, L1.T)[0]).T
        if stf is None:
            b = np.hstack((
                           np.atleast_2d(1),
                           np.zeros((1, order - 1)),
                           np.atleast_2d(1)
                         ))
    if stf is not None:
        # Compute the TF from each feed-in to the output
        # and solve for coefficients which yield the best match
        # THIS CODE IS NOT OPTIMAL,  in terms of computational efficiency.
        stfList = []
        for i in range(0, order + 1):
            bi = np.zeros(1, order + 1)
            bi[i] = 1
            ABCD = stuffABCD(a, g, bi, c, form)
            if form[2:4] == 'FF':
                ABCD[0, order + 1] = -1
            _, stfListi = calculateTF(ABCD)
            stfList.append(stfListi)
        # Build the matrix equation b A  =  x and solve it.
        A = np.zeros((order + 1, max(zSet.shape)))
        for i in range(0, order + 1):
            A[i, :] = evalTF(stfList[i], zSet)
        x = evalTF(stf, zSet)
        x = x[:].T
        b = np.real(x / A)
    return a.squeeze(), g.squeeze(), b.squeeze(), c.squeeze()

