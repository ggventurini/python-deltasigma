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

from __future__ import division, print_function
import numpy as np
from ._partitionABCD import partitionABCD
from ._utils import carray, diagonal_indices

def stuffABCD(a, g, b, c, form='CRFB'):
    """Calculate the ABCD matrix from the parameters of a modulator topology.

    **Parameters:**

    a : array_like
        Feedback/feedforward coefficients from/to the quantizer. Length :math:`n`.
    g : array_like
        Resonator coefficients. Length :math:`floor(n/2)`.
    b : array_like
        Feed-in coefficients from the modulator input to each integrator. Length :math:`n + 1`.
    c : array_like
        Integrator inter-stage coefficients. Length :math:`n`.
    form : str, optional
        See :func:`realizeNTF` for a list of supported structures.

    **Returns:**

    ABCD : ndarray
        A state-space description of the modulator loop filter.
    
    .. seealso:: :func:`mapABCD`, the inverse function.
    """
    # Code common to all structures.
    # a, g, b, and c are internally row vectors
    a = carray(a).reshape((1, -1))
    g = carray(g).reshape((1, -1))
    b = carray(b).reshape((1, -1))
    c = carray(c).reshape((1, -1))
    order = max(a.shape)
    odd = order % 2
    even = 1 - odd # OMG
    ABCD = np.zeros((order + 1, order + 2))
    if np.isscalar(b) or max(b.shape) == 1:
        b = np.atleast_2d(b)
        b = np.hstack((b, np.zeros((1, order))))

    # mutually exclusive matrix stuffing
    if form == 'CRFB':
        # C=(0 0...c_n)
        # This is done as part of the construction of A, below
        # B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD[:, order] = b
        # B2 = -(a_1 a_2... a_n)
        ABCD[:order, order + 1] = -a
        diagonal = diagonal_indices(ABCD[:order, :order])
        ABCD[diagonal] = np.ones((order,))
        subdiag = [i[even:order:2] for i in diagonal_indices(ABCD, -1)]
        ABCD[subdiag] = c[0, even:order:2]
        supdiag = [i[odd:order:2] for i in diagonal_indices(ABCD[:order, :order], +1)]
        ABCD[supdiag] = -g.reshape((-1,))
        # row numbers of delaying integrators
        dly = np.arange(odd + 1, order, 2)
        ABCD[dly, :] = ABCD[dly, :] + np.dot(np.diag(c[0, dly - 1]), ABCD[dly - 1, :])
    elif form == 'CRFF':
        # B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD[:, order] = b
        # B2 = -(c_1 0... 0)
        ABCD[0, order + 1] = -c[0, 0]
        # number of elements = order
        diagonal = diagonal_indices(ABCD[:order, :order])
        ABCD[diagonal] = np.ones((order, ))
        #subdiag = diagonal[:order - 1:2] + 1
        subdiag = [i[:order - 1:2] for i in diagonal_indices(ABCD, -1)]
        ABCD[subdiag] = c[0, 1:order:2]
        if even:
            # rows to have g*(following row) subtracted.
            multg = np.arange(0, order, 2)
            ABCD[multg, :] = ABCD[multg, :] \
                             - np.dot(np.diag(g.reshape((-1,))), ABCD[multg + 1, :])
        else:
            if order >= 3:
                #supdiag = diagonal[2:order:2] - 1
                supdiag = [i[1:order:2] for i in
                              diagonal_indices(ABCD[:order, :order], +1)]
                ABCD[supdiag] = -g.reshape((-1,))
        # rows to have c*(preceding row) added.
        multc = np.arange(2, order, 2)
        ABCD[multc, :] = ABCD[multc, :] + np.dot(np.diag(c[0, multc]), ABCD[multc - 1, :])
        print(order)
        ABCD[order, :order:2] = a[0, :order:2]
        for i in range(1, order, 2):
            ABCD[order, :] = ABCD[order, :] + a[0, i]*ABCD[i, :]
    elif form == 'CIFB':
        # C=(0 0...c_n)
        # This is done as part of the construction of A, below
        # B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD[:, order] = b
        # B2 = -(a_1 a_2... a_n)
        ABCD[:order, order + 1] = -a
        diagonal = diagonal_indices(ABCD[:order, :order])
        ABCD[diagonal] = np.ones((order,))
        subdiag = diagonal_indices(ABCD, -1)
        ABCD[subdiag] = c.reshape((-1,))
        supdiag = [i[odd:order+odd:2] for i in diagonal_indices(ABCD[:order, :order], +1)]
        ABCD[supdiag] = -g.reshape((-1,))
    elif form == 'CIFF':
        # B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD[:, order] = b
        # B2 = -(c_1 0... 0)
        ABCD[0, order + 1] = -c[0, 0]
        diagonal = diagonal_indices(ABCD[:order, :order])
        ABCD[diagonal] = np.ones((order, ))
        subdiag = diagonal_indices(ABCD[:order, :order], -1)
        ABCD[subdiag] = c[0, 1:]
        # C = (a_1 a_2... a_n)
        ABCD[order, :order] = a[0, :order]
        #supdiag = diagonal[odd + 1:order:2] - 1
        supdiag = [i[odd:order+odd:2] for i in  
                      diagonal_indices(ABCD[:order, :order], +1)]
        ABCD[supdiag] = -g.reshape((-1,))
    elif form == 'CRFBD':
        # C=(0 0...c_n)
        ABCD[order, order - 1] = c[0, order - 1]
        # B1 = (b_1 b_2... b_n), D=(b_n+1 0)
        ABCD[:, order] = b
        # B2 = -(a_1 a_2... a_n)
        ABCD[:order, order + 1] = -a
        diagonal = diagonal_indices(ABCD[:order, :order])
        ABCD[diagonal] = np.ones((order, ))
        # row numbers of delaying integrators
        dly = np.arange(odd, order, 2)
        print(dly, diagonal_indices(ABCD[:order, :order], -1))
        subdiag = [np.atleast_1d(i)[dly] for i in
                      diagonal_indices(ABCD[:order, :order], -1)]
        ABCD[subdiag] = c[0, dly]
        supdiag = [np.atleast_1d(i)[dly] for i in 
                      diagonal_indices(ABCD[:order, :order], +1)]
        supdiag = [np.atleast_1d(i[odd:]) for i in supdiag]
        ABCD[dly, :] = ABCD[dly, :] - np.dot(np.diag(g.reshape((-1))), \
                                             ABCD[dly + 1, :])
        if order > 2:
            coupl = np.arange(even + 1, order, 2)
            ABCD[coupl, :] = ABCD[coupl, :] + np.dot(np.diag(c[0, coupl - 1]), \
                             ABCD[coupl - 1, :])
    elif form == 'CRFFD':
        diagonal = diagonal_indices(ABCD[:order, :order])
        subdiag = diagonal_indices(ABCD[:order, :order], -1)
        supdiag = diagonal_indices(ABCD[:order, :order], +1)
        # B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD[:, order] = b.T
        # B2 = -(c_1 0... 0)
        ABCD[0, order + 1] = -c[0, 0]
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
        ABCD[:, order] = b
        ## B2 = -(c_1 0... 0)
        ABCD[0, order + 1] = -c[0, 0]
        diagonal = diagonal_indices(ABCD[:order, :order])
        ABCD[diagonal] = np.ones((order, ))
        subdiag = diagonal_indices(ABCD[:order, :order], -1)
        ABCD[subdiag] = c[0, 1:]
        # code based on case 'CRFF':
        ## rows to have g*(following row) subtracted.
        multg = np.arange(odd, order - 1, 2)
        ABCD[multg, :] = ABCD[multg, :] - np.dot(np.diag(g.reshape((-1,))), \
                                                 ABCD[multg + 1, :])
        # code copied from case 'CIFF':
        ## C = (a_1 a_2... a_n)
        ABCD[order, :order] = a[0, :order]
    else:
        raise TypeError('Form %s is not yet supported.', form)
    return ABCD

def test_stuffABCD():
	"""Test function for stuffABCD()"""
	# we test for the following values:
	#orders = (2, 3, 4, 5)
	#osr = 32
	#nlev = 2
	#f0s = (0., 0.25)
	#Hinf = 1.5
	#forms = ('CRFB', 'CRFF', 'CIFB', 'CIFF', 'CRFBD', 'Stratos')
	tv = {0.:{'CIFB':{2:{'a':(0.2164, 0.7749),
	                      'g':(0, ),
	                      'b':(0.2164, 0.7749, 1.0000),
	                      'c':(1., 1. ),
	                      'ABCD':((1.0000, 0., 0.2164, -0.2164),
	                              (1.0000, 1.0000, 0.7749, -0.7749),
                                  (0, 1.0000, 1.0000, 0))
	                     },
	                   3:{'a':(0.0444, 0.2843, 0.8025),
	                      'g':(0.0058, ),
	                      'b':(0.0444, 0.2843, 0.8025, 1.),
	                      'c':(1., 1., 1.),
	                      'ABCD':((1.0000, 0, 0, 0.0444, -0.0444),
                                  (1.0000, 1.0000, -0.0058, 0.2843, -0.2843),
                                  (0, 1.0000, 1.0000, 0.8025, -0.8025),
                                  (0, 0, 1.0000, 1.0000, 0))
	                     },
	                   4:{'a':(0.0062, 0.0655, 0.3042, 0.8089),
	                      'g':(0., 0.0069),
	                      'b':(0.0062, 0.0655, 0.3042, 0.8089, 1.),
	                      'c':(1., 1., 1., 1),
	                      'ABCD':((1.0000, 0, 0, 0, 0.0062, -0.0062),
                                  (1.0000, 1.0000, 0, 0, 0.0655, -0.0655),
                                  (0, 1.0000, 1.0000, -0.0069, 0.3042, -0.3042),
                                  (0, 0, 1.0000, 1.0000, 0.8089, -0.8089),
                                  (0, 0, 0, 1.0000, 1.0000, 0))
	                     },
	                   5:{'a':(0.0007, 0.0095, 0.0731, 0.3100, 0.81309),
	                      'g':(0.0028, 0.0079),
	                      'b':(0.0007, 0.0095, 0.0731, 0.3100, 0.8130, 1.),
	                      'c':(1., 1., 1., 1., 1.),
	                      'ABCD':((1., 0, 0, 0, 0, 0.0007, -0.0007),
                                  (1., 1.0000, -0.0028, 0, 0, 0.0095, -0.0095),
                                  (0, 1.0000, 1.0000, 0, 0, 0.0731, -0.0731),
                                  (0, 0, 1.0000, 1.0000, -0.0079, 0.3100, -0.3100),
                                  (0, 0, 0, 1.0000, 1.0000, 0.8130, -0.8130),
                                  (0, 0, 0, 0, 1.0000, 1.0000, 0))
	                     }
	                  },
	           'CRFB':{2:{'a':(0.2164, 0.5585),
	                      'g':(0, ),
	                      'b':(0.2164, 0.5585, 1.0000),
	                      'c':(1., 1. ),
	                      'ABCD':((1.0000, 0., 0.2164, -0.2164),
                                  (1.0000, 1.0000, 0.7749, -0.7749),
                                  (0., 1.0, 1.0, 0))
	                     },
	                   3:{'a':(0.0444, 0.2399, 0.5569),
	                      'g':(0.0058, ),
	                      'b':(0.0444, 0.2399, 0.5569, 1.),
	                      'c':(1., 1., 1.),
	                      'ABCD':((1.0000, 0., 0., 0.0444, -0.0444),
                                  (1.0000, 1.0, -0.0058, 0.2399, -0.2399),
                                  (1.0000, 1.0,  0.9942, 0.7967, -0.7967),
                                  (0.,     0.,   1.0000, 1.0000,  0.))
	                     },
	                   4:{'a':(0.0062, 0.0530, 0.2449, 0.5571),
	                      'g':(0, 0.0069),
	                      'b':(0.0062, 0.0530, 0.2449, 0.5571, 1.),
	                      'c':(1., 1., 1., 1),
	                      'ABCD':((1., 0., 0.,  0., 0.0062, -0.0062),
                                  (1., 1., 0.,  0., 0.0592, -0.0592),
                                  (0., 1., 1., -0.0069, 0.2449, -0.2449),
                                  (0., 1., 1.,  0.9931, 0.8020, -0.8020),
                                  (0., 0., 0.,  1.0000, 1.0000,  0.))
	                     },
	                   5:{'a':(0.0007, 0.0084, 0.0550, 0.2443, 0.5579),
	                      'g':(0.0028, 0.0079),
	                      'b':(0.0007, 0.0084, 0.0550, 0.2443, 0.5579, 1.),
	                      'c':(1., 1., 1., 1., 1.),
	                      'ABCD':((1., 0.,  0.,     0., 0., 0.0007, -0.0007),
                                  (1., 1., -0.0028, 0., 0., 0.0084, -0.0084),
                                  (1., 1.,  0.9972, 0., 0., 0.0633, -0.0633),
                                  (0., 0.,  1.0000, 1., -0.0079, 0.2443, -0.2443),
                                  (0., 0.,  1.0000, 1.,  0.9921, 0.8023, -0.8023),
                                  (0., 0.,  0.,     0.,  1.0000, 1.0000,  0.))
	                     }
	                   },
	           'CRFF':{2:{'a':(0.5585, 0.2164),
	                      'g':(0, ),
	                      'b':(1., 0., 1.),
	                      'c':(1., 1. ),
	                      'ABCD':((1., 0., 1., -1.),
                                  (1., 1., 0.,  0.),
                                  (0.7749, 0.2164, 1., 0.))
	                     },
	                   3:{'a':(0.5569, 0.2399, 0.0412),
	                      'g':(0.0058, ),
	                      'b':(1., 0., 0., 1.),
	                      'c':(1., 1., 1.),
	                      'ABCD':((1., 0.,  0., 1., -1.),
                                  (1., 1.,  -0.0058, 0., 0.),
                                  (1., 1.,   0.9942, 0., 0.),
                                  (0.7967, 0.2399, 0.0398, 1.0000, 0.))
	                     },
	                   4:{'a':(0.5571, 0.2449, 0.0492, 0.0046),
	                      'g':(0, 0.0069),
	                      'b':(1., 0., 0., 0., 1.),
	                      'c':(1., 1., 1., 1),
	                      'ABCD':((1., 0., 0, 0, 1.0000, -1.0000),
                                  (1.0000, 1.0000, 0, 0, 0, 0),
                                  (1.0000, 1.0000, 0.9931, -0.0069, 0, 0),
                                  (0, 0, 1.0000, 1.0000, 0, 0),
                                  (0.8020, 0.2449, 0.0537, 0.0046, 1.0000, 0))
	                     },
	                   5:{'a':(0.5579, 0.2443, 0.0505, 0.0071, 0.0003),
	                      'g':(0.0028, 0.0079),
	                      'b':(1., 0., 0., 0., 0., 1.),
	                      'c':(1., 1., 1., 1., 1.),
	                      'ABCD':((1.0000, 0, 0, 0, 0, 1.0000, -1.0000),
                                  (1.0000, 1.0000, -0.0028, 0, 0, 0, 0),
                                  (1.0000, 1.0000, 0.9972,  0, 0, 0, 0),
                                  (0, 0, 1.0000, 1.0000, -0.0079, 0, 0),
                                  (0, 0, 1.0000, 1.0000,  0.9921, 0, 0),
                                  (0.8023, 0.2443, 0.0570, 0.0071, 0.0002, 1., 0))
	                     }
	                   },
	           'CIFF':{2:{'a':(0.7749, 0.2164),
	                      'g':(0., ),
	                      'b':(1., 0., 1.),
	                      'c':(1., 1. ),
	                      'ABCD':((1.0000, 0, 1.0000, -1.0000),
                                  (1.0000, 1.0000, 0,  0),
                                  (0.7749, 0.2164, 1.0000, 0))
	                     },
	                   3:{'a':(0.8025, 0.2843, 0.0398),
	                      'g':(0.0058, ),
	                      'b':(1., 0., 0., 1.),
	                      'c':(1., 1., 1.),
	                      'ABCD':((1.0000, 0, 0, 1.0000, -1.0000),
                                  (1.0000, 1.0000, -0.0058, 0, 0),
                                  (0, 1.0000, 1.0000, 0, 0),
                                  (0.8025, 0.2843, 0.0398, 1.0000, 0))
	                     },
	                   4:{'a':(0.8089, 0.3042, 0.0599, 0.0041),
	                      'g':(0., 0.0069),
	                      'b':(1., 0., 0., 0., 1.),
	                      'c':(1., 1., 1., 1),
	                      'ABCD':((1.0000, 0, 0, 0, 1.0000, -1.0000),
                                  (1.0000, 1.0000, 0, 0, 0, 0),
                                  (0, 1.0000, 1.0000, -0.0069, 0, 0),
                                  (0, 0, 1.0000, 1.0000, 0, 0),
                                  (0.8089, 0.3042, 0.0599, 0.0041, 1.0000, 0))
	                     },
	                   5:{'a':(0.8130, 0.3100, 0.0667, 0.0080, 0.0001),
	                      'g':(0.0028, 0.0079),
	                      'b':(1., 0., 0., 0., 0., 1.),
	                      'c':(1., 1., 1., 1., 1.),
	                      'ABCD':((1.0000, 0, 0, 0, 0, 1.0000, -1.0000),
                                  (1.0000, 1.0000, -0.0028, 0, 0, 0, 0),
                                  (0, 1.0000, 1.0000, 0, 0, 0, 0),
                                  (0, 0., 1.0000, 1.0000, -0.0079, 0, 0),
                                  (0, 0., 0, 1.0000, 1.0000, 0, 0),
                                  (0.8130, 0.3100, 0.0667, 0.0080, 0.0001,  1., 0))
	                     }
	                  },
	           'CRFBD':{2:{'a':(0.2164, 0.7749),
	                       'g':(0, ),
	                       'b':(0.2164, 0.7749, 1.),
	                       'c':(1., 1. ),
	                       'ABCD':((1.0000, 0, 0.2164, -0.2164),
                                   (1.0000, 1.0000, 0.7749, -0.7749),
                                   (0, 1.0000, 1.0000, 0))
	                      },
	                    3:{'a':(0.0444, 0.2399, 0.7967),
	                       'g':(0.0058, ),
	                       'b':(0.0444, 0.2399, 0.7967, 1.),
	                       'c':(1., 1., 1.),
	                       'ABCD':((1.0000, 0, 0, 0.0444, -0.0444),
                                   (1.0000, 0.9942, -0.0058, 0.2796, -0.2796),
                                   (0, 1.0000, 1.0000, 0.7967, -0.7967),
                                   (0, 0, 1.0000, 1.0000, 0))
	                      },
	                    4:{'a':(0.0062, 0.0592, 0.2449, 0.8020),
	                       'g':(0, 0.0069),
	                       'b':(0.0062, 0.0592, 0.2449, 0.8020, 1.),
	                       'c':(1., 1., 1., 1),
	                       'ABCD':((1.0000, 0, 0, 0, 0.0062, -0.0062),
                                   (1.0000, 1.0000, 0, 0, 0.0592, -0.0592),
                                   (1.0000, 1.0000, 0.9931, -0.0069, 0.2987, -0.2987),
                                   (0, 0, 1.0000, 1.0000, 0.8020, -0.8020),
                                   (0, 0, 0, 1.0000, 1.0000, 0))
	                      },
	                    5:{'a':(0.0007, 0.0084, 0.0633, 0.2443, 0.8023),
	                       'g':(0.0028, 0.0079),
	                       'b':(0.0007, 0.0084, 0.0633, 0.2443, 0.8023, 1.),
	                       'c':(1., 1., 1., 1., 1.),
	                       'ABCD':((1.0000, 0, 0, 0, 0, 0.0007, -0.0007),
                                   (1.0000, 0.9972, -0.0028, 0, 0, 0.0089, -0.0089),
                                   (0, 1.0000, 1.0000, 0, 0, 0.0633, -0.0633),
                                   (0, 1., 1., 0.9921, -0.0079, 0.3013, -0.3013),
                                   (0, 0,  0,  1.0000,  1.0000, 0.8023,  -0.8023),
                                   (0, 0,  0,  0, 1.0000, 1.0000, 0))
	                      }
	                   },
	           #'CRFFD':{2:{'a':(),
	           #            'g':(),
	           #            'b':(),
	           #            'c':()
	           #           },
	           #         3:{'a':(),
	           #            'g':(),
	           #            'b':(),
	           #            'c':()
	           #           },
	           #         4:{'a':(),
	           #            'g':(),
	           #            'b':(),
	           #            'c':()
	           #           },
	           #         5:{'a':(),
	           #            'g':(),
	           #            'b':(),
	           #            'c':()
	           #           }
	           #        }
	           'Stratos':{2:{'a':(0.7749, 0.2164),
	                       'g':(0, ),
	                       'b':(1, 0., 1.),
	                       'c':(1., 1. ),
	                       'ABCD':((1.0000, 0, 1.0000, -1.0000),
                                   (1.0000, 1.0000, 0, 0),
                                   (0.7749, 0.2164, 1.0000, 0))
	                      },
	                    3:{'a':(0.7967, 0.2796, 0.0398),
	                       'g':(0.0058, ),
	                       'b':(1., 0., 0., 1.),
	                       'c':(1., 1., 1.),
	                       'ABCD':((1.0000, 0, 0, 1.0000, -1.0000),
                                   (1.0000, 0.9942, -0.0058, 0, 0),
                                   (0, 1.0000, 1.0000, 0, 0),
                                   (0.7967, 0.2796, 0.0398, 1.0000, 0))
	                      },
	                    4:{'a':(0.8020, 0.2987, 0.0579, 0.0042),
	                       'g':(0, 0.0069),
	                       'b':(1., 0., 0., 0., 1.),
	                       'c':(1., 1., 1., 1),
	                       'ABCD':((1.0000, 0, 0, 0, 1.0000, -1.0000),
                                   (1.0000, 1.0000, 0, 0, 0, 0),
                                   (0, 1.0000, 0.9931, -0.0069, 0, 0),
                                   (0, 0, 1.0000, 1.0000, 0, 0),
                                   (0.8020, 0.2987, 0.0579, 0.0042, 1.0000, 0))
	                      },
	                    5:{'a':(0.8023, 0.3013, 0.0643, 0.0075, 0.0001),
	                       'g':(0.0028, 0.0079),
	                       'b':(1., 0., 0., 0., 0., 1.),
	                       'c':(1., 1., 1., 1., 1.),
	                       'ABCD':((1., 0, 0, 0, 0, 1.0000, -1.0000),
                                   (1.0000, 0.9972, -0.0028, 0, 0, 0, 0),
                                   (0, 1.0000, 1.0000, 0, 0, 0, 0),
                                   (0, 0, 1.0000, 0.9921, -0.0079, 0, 0),
                                   (0, 0, 0, 1.0000, 1.0000, 0, 0),
                                   (0.8023, 0.3013, 0.0643, 0.0075, 0.0001, 1., 0))
	                      }
	                   }
	          },
	      0.25:{'CIFB':{2:{'a':(0.3333, 2.0000),
	                      'g':(1., ),
	                      'b':(0.3333, 2., 1.0000),
	                      'c':(1., 1. ),
	                      'ABCD':((1.0000, -1.0000, 0.3333, -0.3333),
                                  (1.0000,  1.0000, 2.0000, -2.0000),
                                  (0.,      1.0000, 1.0000,  0.))
	                     },
	                   4:{'a':(-3.5585, 2.4503, 5.22512, 4.0000),
	                      'g':(1., 1.),
	                      'b':(-3.5585, 2.4503, 5.22512, 4.0000, 1.),
	                      'c':(1., 1., 1., 1),
	                      'ABCD':((1.0000, -1.0000, 0, 0, -3.5585,  3.5585),
                                  (1.0000,  1.0000, 0, 0,  2.4503, -2.4503),
                                  (0, 1.0000, 1.0000, -1.0000, 5.2251, -5.2251),
                                  (0, 0, 1.0000, 1.0000, 4.0000, -4.0000),
                                  (0, 0, 0, 1.0000, 1.0000, 0))
	                     }
	                  },
	           'CRFB':{2:{'a':(-0.6667, 0.6667),
	                      'g':(2.0, ),
	                      'b':(-0.6667, 0.6667, 1.0000),
	                      'c':(1., 1. ),
	                      'ABCD':((1.0000, -2.0000, -0.6667, 0.6667),
                                  (1.0000, -1.0000, -0.0000, 0.0000),
                                   (0., 1.0000, 1.0000, 0.))
	                     },
	                   4:{'a':(-0.2164, 0., -0.5585, 0.5585),
	                      'g':(2.0, 2.0),
	                      'b':(-0.2164, 0.0, -0.5585, 0.5585, 1.),
	                      'c':(1., 1., 1., 1),
	                      'ABCD':((1., -2., 0, 0, -0.2164, 0.2164),
                                  (1., -1., 0, 0, -0.2164, 0.2164),
                                  (0., 1., 1., -2., -0.5585, 0.5585),
                                  (0, 1., 1., -1., 0., 0.),
                                  (0., 0., 0., 1., 1., 0.))
	                     }
	                   },
	           'CRFF':{2:{'a':(0.6667, -0.6667),
	                      'g':(2.0, ),
	                      'b':(1., 0., 1.),
	                      'c':(1., 1. ),
	                      'ABCD':((-1., -2., 1., -1.),
                                  (1., 1., 0., 0.),
                                  (0., -0.6667, 1., 0))
	                     },
	                   4:{'a':(0.5585, -0.5585, 0., -0.2164),
	                      'g':(2., 2.),
	                      'b':(1., 0., 0., 0., 1.),
	                      'c':(1., 1., 1., 1),
	                      'ABCD':((-1., -2., 0., 0., 1., -1.),
                                  ( 1., 1., 0., 0., 0., 0.),
                                  (1., 1., -1., -2., 0., 0.),
                                  (0., 0., 1., 1., 0., 0.),
                                  (0., -0.5585, -0.2164, -0.2164, 1., 0.))
	                     }
	                   },
	           'CIFF':{2:{'a':(2., 0.3333),
	                      'g':(1., ),
	                      'b':(1., 0., 1.),
	                      'c':(1., 1. ),
	                      'ABCD':((1., -1., 1., -1.),
                                  (1., 1., 0., 0.),
                                  (2., 0.3333, 1., 0.))
	                     },
	                   4:{'a':(4., 5.2251, 2.4503, -3.5585),
	                      'g':(1., 1.),
	                      'b':(1., 0., 0., 0., 1.),
	                      'c':(1., 1., 1., 1),
	                      'ABCD':((1., -1., 0., 0., 1., -1.),
                                  (1.,  1., 0., 0., 0., 0.),
                                  (0.,  1., 1., -1., 0., 0),
                                  (0., 0., 1., 1., 0., 0),
                                  (4., 5.2251, 2.4503, -3.5585, 1., 0.))
	                     }
	                  },
	           'CRFBD':{2:{'a':(-0.6667, 0.),
	                       'g':(2.0, ),
	                       'b':(-0.6667, 0., 1.),
	                       'c':(1., 1. ),
	                       'ABCD':((-1.0000, -2.0000, -0.6667, 0.6667),
                                   (1.0000, 1.0000, -0.0000, 0.0000),
                                   (0., 1.0000, 1.0000, 0.))
	                      },
	                    4:{'a':(-0.2164, -0.2164, -0.5585, 0.),
	                       'g':(2., 2.),
	                       'b':(-0.2164, -0.2164, -0.5585, 0., 1.),
	                       'c':(1., 1., 1., 1),
	                       'ABCD':((-1.0000, -2.0000, 0, 0, 0.2164, -0.2164),
                                   (1.0000, 1.0000, 0, 0, -0.2164, 0.2164),
                                   (1., 1., -1., -2.0000, -0.7749, 0.7749),
                                   (0., 0., 1., 1., 0., 0.),
                                   (0., 0., 0., 1., 1., 0.))
	                      }
	                   },
	           #'CRFFD':{2:{'a':(),
	           #            'g':(),
	           #            'b':(),
	           #            'c':()
	           #           },
	           #         4:{'a':(),
	           #            'g':(),
	           #            'b':(),
	           #            'c':()
	           #           }
	           #        }
	           'Stratos':{2:{'a':(0., -0.6667),
	                       'g':(2., ),
	                       'b':(1, 0., 1.),
	                       'c':(1., 1. ),
	                       'ABCD':((-1., -2., 1., -1.),
                                   (1., 1., 0., 0.),
                                   (-0., -0.6667, 1., 0.))
	                      },
	                    4:{'a':(0., -0.7749, 0., 0.2164),
	                       'g':(2.0, 2.0),
	                       'b':(1., 0., 0., 0., 1.),
	                       'c':(1., 1., 1., 1),
	                       'ABCD':((-1., -2., 0., 0., 1., -1.),
                                   (1., 1., 0., 0., 0., 0.),
                                   (0., 1., -1., -2., 0., 0.),
                                   (0., 0., 1., 1., 0., 0.),
                                   (0., -0.7749, 0., 0.2164, 1., 0.))
	                      }
	                   }
	          }
	      }

	for f0 in tv:
		for form in tv[f0]:
			for order in tv[f0][form]:
				# Optimized zero placement
				print("Testing form: %s, order: %d, f0: %f" % \
				      (form, order, f0))
				a = np.array(tv[f0][form][order]['a']).reshape((1, -1))
				g = np.array(tv[f0][form][order]['g']).reshape((1, -1))
				b = np.array(tv[f0][form][order]['b']).reshape((1, -1))
				c = np.array(tv[f0][form][order]['c']).reshape((1, -1))
				ABCD = stuffABCD(a, g, b, c, form)
				print(ABCD)
				assert np.allclose(ABCD, tv[f0][form][order]['ABCD'], 
				            atol=1e-4, rtol=1e-3)
	return 
