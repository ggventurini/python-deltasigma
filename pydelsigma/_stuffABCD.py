# -*- coding: utf-8 -*-
# _stuffABCD.py
# Module providing the stuffABCD function
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

"""Module providing the stuffABCD() function
"""

from __future__ import division
import numpy as np
from ._partitionABCD import partitionABCD

def stuffABCD(a, g, b, c, form='CRFB'):
    """Compute the ABCD matrix for the specified structure.
    See realizeNTF.m for a list of supported structures.
    
    mapABCD is the inverse function.
    """
    # Code common to all structures.
    order = max(a.shape)
    odd = order % 2
    even = 1 - odd # OMG
    ABCD = np.zeros((order + 1, order + 2))
    if np.isscalar(b) or max(b.shape) == 1:
        b = np.atleast_2d(b)
        b = np.vstack((b, np.zeros(1, order)))

    # mutually exclusive matrix stuffing
    if 'CRFB' == form:
        # C=(0 0...c_n)
        # This is done as part of the construction of A, below
        # B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD[:, order] = b.T
        # B2 = -(a_1 a_2... a_n)
        ABCD[:order, order + 1] = -a.T
        # diagonal in 2d:
        # :-1, 0 is the 2d version of 1:(order+2):order*(order+1) 
        # subdiag in 2d
        # even:-1:2, 1 is the 2d version of diagonal((1+even):2:order)+1
        # supdiag in 2d
        # odd:-1:2, -1 is the 2d version of subdiag((1+odd):length(subdiag))-2
        ABCD[:-1, 0] = np.ones((1, order))
        ABCD[even:-1:2, 1] = c[0, even:order:2]
        #print range(odd, order
        ABCD[odd:-1:2, -1] = -g
        # row numbers of delaying integrators
        dly = np.arange(odd + 1, order, 2)
        ABCD[dly, :] = ABCD[dly, :] + np.dot(np.diag(c[0, dly - 1]), ABCD[dly - 1, :])
    elif 'CRFF' == form:
        # B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD[:, order] = b.T
        # B2 = -(c_1 0... 0)
        ABCD[0, order + 1] = -c[0]
        # number of elements = order
        diagonal = np.arange(0, order*(order + 1), order + 2)
        ABCD[diagonal] = np.ones((1, order))
        subdiag = diagonal[:order - 1:2] + 1
        ABCD[subdiag] = c[1:order:2]
        if even:
            # rows to have g*(following row) subtracted.
            multg = range(0, order, 2)
            ABCD[multg, :] = ABCD[multg, :] - np.diag(g) * ABCD[multg + 1, :]
        else:
            if order >= 3:
                supdiag = diagonal[2:order:2] - 1
                ABCD[supdiag] = -g
        # rows to have c*(preceding row) added.
        multc = range(2, order, 2)
        ABCD[multc, :] = ABCD[multc, :] + np.diag(c[multc]) * ABCD[multc - 1, :]
        ABCD[order, :order:2] = a[:order:2]
        for i in range(1, order, 2):
            ABCD[order, :] = ABCD[order, :] + a[i]*ABCD[i, :]
    elif 'CIFB' == form:
        # C=(0 0...c_n)
        # This is done as part of the construction of A, below
        # B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD[:, order] = b.T
        # B2 = -(a_1 a_2... a_n)
        ABCD[:order, order + 1] = -a.T
        diagonal = np.arange(0, order*(order + 1), order + 2)
        ABCD[diagonal] = np.ones((1, order))
        subdiag = diagonal[:order] + 1
        ABCD[subdiag] = c
        supdiag = diagonal[odd + 1:order+odd:2] - 1
        ABCD[supdiag] = -g
    elif 'CIFF' == form:
        # B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD[:, order] = b.T
        # B2 = -(c_1 0... 0)
        ABCD[0, order + 1] = -c[0]
        diagonal = np.arange(0, order*(order + 1), order + 2)
        ABCD[diagonal] = np.ones((1, order))
        subdiag = diagonal[:order - 1] + 1
        ABCD[subdiag] = c[1:]
        # C = (a_1 a_2... a_n)
        ABCD[order, :order] = a[:order]
        supdiag = diagonal[odd + 1:order:2] - 1
        ABCD[supdiag] = -g
    elif 'CRFBD' == form:
        # C=(0 0...c_n)
        ABCD[order, order - 1] = c[order - 1]
        # B1 = (b_1 b_2... b_n), D=(b_n+1 0)
        ABCD[:, order] = b.T
        # B2 = -(a_1 a_2... a_n)
        ABCD[:order, order + 1] = -a.T
        diagonal = np.arange(0, order*(order + 1), order + 2)
        ABCD[diagonal] = np.ones((1, order))
        # row numbers of delaying integrators
        dly = np.arange(odd, order, 2)
        subdiag = diagonal[dly] + 1
        ABCD[subdiag] = c[dly]
        #supdiag = subdiag[odd - 1:max(subdiag.shape)] - 2
        ABCD[dly, :] = ABCD[dly, :] - np.diag(g)*ABCD[dly + 1, :]
        if order > 2:
            coupl = np.arange(even + 1, order, 2)
            ABCD[coupl, :] = ABCD[coupl, :] + np.diag(c[coupl - 1]) * \
                             ABCD[coupl - 1, :]
    elif 'CRFFD' == form:
        diagonal = np.arange(0, order*(order + 1), order + 2)
        subdiag = diagonal[:order - 1] + 1
        supdiag = diagonal[1:order] - 1
        # B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD[:, order] = b.T
        # B2 = -(c_1 0... 0)
        ABCD[0, order + 1] = -c[0]
        ABCD[diagonal] = np.ones((1, order))
        # rows to have c*(preceding row) added.
        multc = np.arange(1, order, 2)
        if order > 2:
            ABCD[subdiag[1::2]] = c[2::2]
        if even:
            ABCD[supdiag[::2]] = -g
        else:
	        # subtract g*(following row) from the multc rows
            ABCD[multc, :] = ABCD[multc, :] - np.diag(g) * ABCD[multc + 1, :]
        ABCD[multc, :] = ABCD[multc, :] + np.diag(c[multc]) * ABCD[multc - 1, :]
        # C
        ABCD[order, 1:order:2] = a[1:order:2]
        for i in range(0, order, 2):
            ABCD[order, :] = ABCD[order, :] + a[i] * ABCD[i, :]
        # The above gives y(n+1); need to add a delay to get y(n).
        # Do this by augmenting the states. Note: this means that
        # the apparent order of the NTF is one higher than it actually is.
        A, B, C, D = partitionABCD(ABCD, 2)
        A = np.array([A, np.zeros((order, 1)), C, 0]).reshape(1, -1)
        B = np.vstack((B, D))
        C = np.hstack((np.zeros((1, order)), np.ones((1, 1))))
        D = np.array([0, 0]).reshape(1, 2)
        ABCD = np.array([A, B, C, D]).reshape(1, -1)
    elif 'PFF' == form:
        #B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD[:, order] = b.T
        odd_1 = odd     # !! Bold assumption !!
        odd_2 = 0.      # !! Bold assumption !! 
        gc = g.dot(c[odd_1::2])
        theta = np.arccos(1 - gc/2.)
        if odd_1:
            theta0 = 0.
        else:
            theta0 = theta[0]
        order_2 = 2*max(np.flatnonzero(np.abs(theta - theta0) > 0.5).shape)
        order_1 = order - order_2
        # B2 = -(c_1 0...0 c_n 0...0)
        ABCD[0, order + 1] = -c[0]
        ABCD[order_1, order + 1] = -c[order_1]
        # number of elements = order
        diagonal = np.arange(0, order*(order + 1), order + 2)
        ABCD[diagonal] = np.ones((1, order))
        i = np.hstack((np.arange(0, order_1 - 1, 2),
                       np.arange(order - order_2, order, 2)))
        subdiag = diagonal[i] + 1
        ABCD[subdiag] = c[i + 1]
        if odd_1:
            if order_1 >= 3:
                supdiag = diagonal[2:order_1:2] - 1
                ABCD[supdiag] = -g[:(order_1 - 1)/2.]
        else:
            # rows to have g*(following row) subtracted.
            multg = np.arange(0, order_1, 2)
            ABCD[multg, :] = ABCD[multg, :] - np.diag(g[:order_1/2.]) * \
                             ABCD[multg + 1, :]
        if odd_2:
            if order_2 >= 3:
                supdiag = diagonal[order_1 + 1:order:2] - 1
                ABCD[supdiag] = -g[:(order_1 - 1)/2.]
        else:
            # rows to have g*(following row) subtracted.
            multg = np.arange(order_1, order, 2)
            gg = g[(order_1 - odd_1)/2:]
            ABCD[multg, :] = ABCD[multg, :] - np.diag(gg)*ABCD[multg + 1, :]
        # Rows to have c*(preceding row) added.
        multc = np.hstack((np.arange(2, order_1, 2), 
                           np.arange(order_1 + 2, order, 2)
                         )).reshape(1, -1)
        ABCD[multc, :] = ABCD[multc, :] + np.diag(c[multc])*ABCD[multc - 1, :]
        # C portion of ABCD
        i = np.hstack((
                       np.arange(0, order_1, 2),
                       np.arange(order_1, order, 2)
                     )).reshape(1,-1)
        ABCD[order, i] = a[i]
        for i in np.hstack((np.arange(1, order_1, 2),
                            np.arange(order_1 + 1, order, 2)
                          )):
            ABCD[order, :] = ABCD[order, :] + a[i] * ABCD[i, :]
    elif 'Stratos' == form:
        # code copied from case 'CIFF':
        ## B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD[:, order] = b.T
        ## B2 = -(c_1 0... 0)
        ABCD[0, order + 1] = -c[0]
        diagonal = np.arange(0, order*(order + 1), order + 2)
        ABCD[diagonal] = np.ones((1, order))
        subdiag = diagonal[:order - 1] + 1
        ABCD[subdiag] = c[1:]
        # code based on case 'CRFF':
        ## rows to have g*(following row) subtracted.
        multg = np.arange(odd, order - 1, 2)
        ABCD[multg, :] = ABCD[multg, :] - np.diag(g) * ABCD[multg + 1, :]
        # code copied from case 'CIFF':
        ## C = (a_1 a_2... a_n)
        ABCD[order, :order] = a[:order]
    else:
        raise TypeError('Form %s is not yet supported.', form)
    return ABCD
