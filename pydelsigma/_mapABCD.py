# -*- coding: utf-8 -*-
# _mapABCD.py
# Module providing the mapABCD function
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

"""Module providing the mapABCD() function
"""

from __future__ import division
import numpy as np

from ._utils import diagonal_indices

def mapABCD(ABCD, form='CRFB'):
    """Compute the coefficients for the specified structure.
    See realizeNTF.m for a list of supported structures.
    stuffABCD is the inverse function.
    """
    order = ABCD.shape[0] - 1
    odd = order % 2
    even = 1 - odd
    diagonal = diagonal_indices(ABCD)
    subdiag = diagonal_indices(ABCD, -1)
    supdiag = map(lambda a: a[odd:order - 1:2], diagonal_indices(ABCD, +1))
    if form in ('CRFB', 'CIFB', 'CRFBD'):
        c = ABCD[subdiag]
        g = -ABCD[supdiag]
        if form == 'CRFB':
            dly = np.arange(1 + odd, order, 2)
            ABCD[dly, :] = ABCD[dly, :] \
                           - np.dot(np.diag(c[dly - 1]), ABCD[dly - 1, :])
        elif form == 'CRFBD':
            dly = np.arange(odd, order, 2)
            ABCD[dly, :] = ABCD[dly, :] + np.dot(np.diag(g), ABCD[dly + 1, :])
            if order > 2:
                coupl = np.arange(1 + even, order, 2)
                ABCD[coupl, :] = ABCD[coupl, :] \
                                 - np.dot(np.diag(c[coupl - 1]), ABCD[coupl - 1, :])
        a = -ABCD[:order, order + 1].T
        b = ABCD[:, order].T
    elif form == 'CRFF':
        a = np.zeros((order,))
        c = np.concatenate((
                            np.array((-ABCD[0, order + 1],)), 
                            ABCD[subdiag][:-1]
                          ))
        g = -ABCD[supdiag]
        if even:
            multg = np.arange(0, order, 2)
            ABCD[multg, :] = ABCD[multg, :] + np.dot(np.diag(g), ABCD[multg + 1, :])
        multc = np.arange(2, order, 2)
        ABCD[multc, :] = ABCD[multc, :] - np.dot(np.diag(c[multc]), ABCD[multc - 1, :])
        a[1:order:2] = ABCD[order, 1:order:2]
        for i in range(1, order, 2):
            ABCD[order, :] = ABCD[order, :] - a[i]*ABCD[i, :]
        a[:order:2] = ABCD[order, :order:2]
        b = ABCD[:, order].T
    elif form == 'CRFFD':
        #order=order - 1
        #odd=rem(order,2)
        #even=not  odd
        diagonal = diagonal_indices(ABCD[:order, :order])
        subdiag = diagonal_indices(ABCD[:order, :order], -1)
        supdiag = map(lambda a: a[1+odd:order:2], 
                      diagonal_indices(ABCD[:order, :order], -1))
        g = -ABCD[supdiag]
        c = np.vstack((-ABCD[0, order + 2], ABCD[subdiag]))
        a = np.zeros((1, order))
        for i in range(0, order, 2):
            a[i] = ABCD[order, i]
            ABCD[order, :] = ABCD[order, :] - a[i]*ABCD[i, :]
        a[1:order:2] = ABCD[order, 1:order:2]
        b = ABCD[:order + 1, order + 1].T
        for i in range(1, order, 2):
            b[i] = b[i] - c[i]*b[i - 1]
            if odd:
                b[i] = b[i] + g[(i - 1)/2]*b[i + 1]
        yscale = ABCD[order + 1, order]
        a = a*yscale
        b[-1] = b[-1]*yscale
    elif form in ('CIFF', 'Stratos'):
        a = ABCD[order, :order]
        c = np.vstack([-ABCD[0, order + 1], ABCD[subdiag[:-1]]])
        g = -ABCD[supdiag]
        b = ABCD[:, order].T
    else:
        raise ValueError('Form %s is not yet supported.' % form)

    return a.squeeze(), g.squeeze(), b.squeeze(), c.squeeze()
