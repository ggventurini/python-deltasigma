# -*- coding: utf-8 -*-
# _realizeNTF_ct.py
# Module providing the realizeNTF_ct function
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

"""Module providing the realizeNTF_ct() function
"""

from __future__ import division, print_function
from warnings import warn

import numpy as np
import numpy.linalg as linalg
from scipy.signal import ss2zpk, dimpulse

from ._evalTFP import evalTFP
from ._impL1 import impL1
from ._padb import padb
from ._pulse import pulse
from ._utils import carray, eps, _get_zpk

def realizeNTF_ct(ntf, form='FB', tdac=(0, 1), ordering=None, bp=None,
                  ABCDc=None, method='LOOP'):
    """Realize an NTF with a continuous-time loop filter.

    **Parameters:**

    ntf : object
        A noise transfer function (NTF).

    form : str, optional
        A string specifying the topology of the loop filter.

         * 'FB': Feedback form,
         * 'FF': Feedforward form

        For the FB structure, the elements of ``Bc`` are calculated
        so that the sampled pulse response matches the L1 impulse
        response.  For the FF structure, ``Cc`` is calculated.

    tdac : sequence, optional
        The timing for the feedback DAC(s). If ``tdac[0] >= 1``,
        direct feedback terms are added to the quantizer.

        Multiple timings (one or more per integrator) for the FB
        topology can be specified by making tdac a list of lists,
        e.g. ``tdac = [[1, 2], [1, 2], [[0.5, 1], [1, 1.5]], []]``

        In this example, the first two integrators have
        DACs with ``[1, 2]`` timing, the third has a pair of
        DACs, one with ``[0.5, 1]`` timing and the other with
        ``[1, 1.5]`` timing, and there is no direct feedback
        DAC to the quantizer.

    ordering : sequence, optional
        A vector specifying which NTF zero-pair to use in each resonator
        Default is for the zero-pairs to be used in the order specified
        in the NTF.

    bp : sequence, optional
        A vector specifying which resonator sections are bandpass.
        The default (``zeros(...)``) is for all sections to be lowpass.

    ABCDc : ndarray, optional
        The loop filter structure, in state-space form.
        If this argument is omitted, ABCDc is constructed according
        to "form."

    method : str, optional
        The default fitting method is ``'LOOP'``, which means that
        the DT and CT loop responses will be matched.
        Alternatively, it is possible to set the method to ``'NTF'``,
        which will result in the NTF responses to be matched.
        See :ref:`discrete-time-to-continuous-time-mapping` for a
        more in-depth discussion.

    **Returns:**

    ABCDc : ndarray
        A state-space description of the CT loop filter

    tdac2 : ndarray
        A matrix with the DAC timings, including ones
        that were automatically added.

    **Example:**

    Realize the NTF :math:`(1 - z^{-1})^2` with a CT system (cf with the
    example at :func:`mapCtoD`).::

        from deltasigma import *
        ntf = ([1, 1], [0, 0], 1)
        ABCDc, tdac2 = realizeNTF_ct(ntf, 'FB')

    Returns:

    ABCDc::

        [[ 0.          0.          1.         -1.        ]
         [ 1.          0.          0.         -1.49999999]
         [ 0.          1.          0.          0.        ]]

    tdac2::

        [[-1. -1.]
         [ 0.  1.]]

    """
    ntf_z, ntf_p, _ = _get_zpk(ntf)
    ntf_z = carray(ntf_z)
    ntf_p = carray(ntf_p)
    order = max(ntf_p.shape)
    order2 = int(np.floor(order/2.))
    odd = order - 2*order2
    # compensate for limited accuracy of zero calculation
    ntf_z[np.abs(ntf_z - 1) < eps**(1./(1. + order))] = 1.
    method = method.upper()
    if method not in ('LOOP', 'NTF'):
        raise ValueError('Unimplemented matching method %s.' % method)
    # check if multiple timings mode
    if (type(tdac) == list or type(tdac) == tuple) and len(tdac) and \
       (type(tdac[0]) == list or type(tdac[0]) == tuple):
        if len(tdac) != order + 1:
            msg = 'For multi-timing tdac, len(tdac) ' + \
                  ' must be order+1.'
            raise ValueError(msg)
        if form != 'FB':
            msg = "Currently only supporting form='FB' " + \
                  'for multi-timing tdac'
            raise ValueError(msg)
        multi_timing = True
    else: # single timing
        tdac = carray(tdac)
        if np.prod(tdac.shape) != 2:
            msg = 'For single-timing tdac, len(tdac) must be 2.'
            raise ValueError(msg)
        tdac.reshape((2,))
        multi_timing = False
    if ordering is None:
        ordering = np.arange(order2)
    if bp is None:
        bp = np.zeros((order2,))
    if not multi_timing:
        # Need direct terms for every interval of memory in the DAC
        n_direct = np.ceil(tdac[1]) - 1
        if tdac[0] > 0 and tdac[0] < 1 and tdac[1] > 1 and tdac[1] < 2:
            n_extra = n_direct - 1 #  tdac pulse spans a sample point
        else:
            n_extra = n_direct
        tdac2 = np.vstack(
                 (np.array((-1, -1)),
                  np.array(tdac).reshape((1, 2)),
                  0.5*np.dot(np.ones((n_extra, 1)), np.array([[-1, 1]]))
                  + np.cumsum(np.ones((n_extra, 2)), 0) + (n_direct - n_extra)
                 ))
    else:
        n_direct = 0
        n_extra = 0
    if ABCDc is None:
        ABCDc = np.zeros((order + 1, order + 2))
        # Stuff the A portion
        if odd:
            ABCDc[0, 0] = np.real(np.log(ntf_z[0]))
            ABCDc[1, 0] = 1
        dline = np.array([0, 1, 2])
        for i in range(order2):
            n = bp[i]
            i1 = 2*i + odd
            zi = 2*ordering[i] + odd
            w = np.abs(np.angle(ntf_z[zi]))
            ABCDc[i1 + dline, i1] = np.array([0, 1, n])
            ABCDc[i1 + dline, i1 + 1] = np.array([-w**2, 0, 1 - n])
        ABCDc[0, order] = 1
        # 2006.10.02 Changed to -1 to make FF STF have +ve gain at DC
        ABCDc[0, order + 1] = -1
    Ac = ABCDc[:order, :order]
    if form == 'FB':
        Cc = ABCDc[order, :order].reshape((1, -1))
        if not multi_timing:
            Bc = np.hstack((np.eye(order), np.zeros((order, 1))))
            Dc = np.hstack((np.zeros((1, order)), np.array([[1]])))
            tp = np.tile(np.array(tdac).reshape((1, 2)), (order + 1, 1))
        else: #Assemble tdac2, Bc and Dc
            tdac2 = np.array([[-1, -1]])
            Bc = None
            Dc = None
            Bci = np.hstack((np.eye(order), np.zeros((order, 1))))
            Dci = np.hstack((np.zeros((1, order)), np.array([[1]])))
            for i in range(len(tdac)):
                tdi = tdac[i]
                if (type(tdi) in (tuple, list)) and len(tdi) and \
                   (type(tdi[0]) in (list, tuple)):
                    for j in range(len(tdi)):
                        tdj = tdi[j]
                        tdac2 = np.vstack((tdac2,
                                           np.array(tdj).reshape(1,-1)))
                        if Bc is not None:
                            Bc = np.hstack((Bc, Bci[:, i].reshape((-1, 1))))
                        else:
                            Bc = Bci[:, i].reshape((-1, 1))
                        if Dc is not None:
                            Dc = np.hstack((Dc, Dci[:, i].reshape((-1, 1))))
                        else:
                            Dc = Dci[:, i].reshape((-1, 1))
                elif len(tdi): # we got tdac[i] = [a, b] where a, b are scalars
                    tdac2 = np.vstack((tdac2,
                                       np.array(tdi).reshape(1,-1)))
                    if Bc is not None:
                        Bc = np.hstack((Bc, Bci[:, i].reshape((-1, 1))))
                    else:
                        Bc = Bci[:, i].reshape((-1, 1))
                    if Dc is not None:
                        Dc = np.hstack((Dc, Dci[:, i].reshape((-1, 1))))
                    else:
                        Dc = Dci[:, i].reshape((-1, 1))
            tp = tdac2[1:, :]
    elif form == 'FF':
        Cc = np.vstack((np.eye(order), np.zeros((1, order))))
        Bc = np.vstack((np.array([[-1]]), np.zeros((order-1, 1))))
        Dc = np.vstack((np.zeros((order, 1)), np.array([[1]])))
        tp = tdac #  2008-03-24 fix from Ayman Shabra
    else:
        raise ValueError('Sorry, no code for form "%s".', form)

    n_imp = np.ceil(2*order + np.max(tdac2[:, 1]) + 1)
    if method == 'LOOP':
        # Sample the L1 impulse response
        y = impL1(ntf, n_imp)
    else:
        # Sample the NTF impulse response
        y = dimpulse((ntf_z, ntf_p, 1., 1.), t=np.arange(n_imp+1))[1][0]
        y = np.atleast_1d(y.squeeze())
    sys_c = []
    for i in range(Bc.shape[1]): # number of inputs
        sys_tmp = []
        for j in range(Cc.shape[0]): # number of outputs
            sys_tmp.append(ss2zpk(Ac, Bc, Cc[j, :], Dc[j, :], input=i))
        sys_c.append(sys_tmp)
    yy = pulse(sys_c, tp, 1, n_imp, 1)
    yy = np.squeeze(yy)
    # Endow yy with n_extra extra impulses.
    # These will need to be implemented with n_extra extra DACs.
    # !! Note: if t1=int, matlab says pulse(sys) @t1 ~=0
    # !! This code corrects this problem.
    if n_extra > 0:
        y_right = padb(np.vstack((np.zeros((1, n_direct)),
                                  np.eye(n_direct))),
                       n_imp + 1)
        # Replace the last column in yy with an ordered set of impulses
        if (n_direct > n_extra):
            yy = np.hstack((yy, y_right[:, 1:]))
        else:
            yy = np.hstack((yy[:, :-1], y_right))

    if method == 'NTF':
        # convolve CT loop response and NTF response
        yynew = None
        for i in range(yy.shape[1]):
            yytmp = np.convolve(yy[:, i], y)[:-n_imp]
            if yynew is None:
                yynew = yytmp.reshape((-1, 1))
            else:
                yynew = np.hstack((yynew, yytmp.reshape((-1, 1))))
        yy = yynew
        e1 = np.zeros(y.shape)
        e1[0] = 1.
        y = y - e1
    # Solve for the coefficients
    x = linalg.lstsq(yy, y)[0]
    if linalg.norm(np.dot(yy, x) - y) > 0.0001:
        warn('Pulse response fit is poor.')
    if form == 'FB':
        if not multi_timing:
            Bc2 = np.hstack((x[:order].reshape((-1, 1)),
                             np.zeros((order, n_extra))))
            if n_extra > 0:
                Dc2 = np.hstack((np.array([[0]]),
                                 x[order:].reshape((-1, 1))))
            else:
                Dc2 = x[order:].reshape((-1, 1))
        else:
            BcDc = np.vstack((Bc, Dc))
            i = np.nonzero(BcDc)
            BcDc[i] = x
            Bc2 = BcDc[:-1, :]
            Dc2 = BcDc[-1, :]
    elif form == 'FF':
        Bc2 = np.hstack((Bc, np.zeros((order, n_extra))))
        Cc = x[:order].reshape((1, -1))
        if n_extra > 0:
            Dc2 = np.hstack((np.array([[0]]), x[order:].T))
        else:
            Dc2 = x[order:].T

    Dc1 = np.zeros((1, 1))
    Dc = np.hstack((Dc1, np.atleast_2d(Dc2)))
    Bc1 = np.vstack((np.ones((1, 1)), np.zeros((order - 1, 1))))
    Bc = np.hstack((Bc1, Bc2))
    # Scale Bc1 for unity STF magnitude at f0
    fz = np.angle(ntf_z)/(2*np.pi)
    f1 = fz[0]
    ibz = np.abs(fz - f1) <= np.abs(fz + f1)
    fz = fz[ibz]
    f0 = np.mean(fz)
    if np.min(np.abs(fz)) < 3*np.min(np.abs(fz - f0)):
        f0 = 0
    L0c = ss2zpk(Ac, Bc1, Cc, Dc1)
    G0 = evalTFP(L0c, ntf, f0)
    if f0 == 0:
        Bc[:, 0] = np.dot(Bc[:, 0],
                          np.abs(np.dot(Bc[0, 1:],
                                        (tdac2[1:, 1] - tdac2[1:, 0]))
                                 /Bc[0, 0]))
    else:
        Bc[:, 0] = Bc[:, 0]/np.abs(G0)
    ABCDc = np.vstack((
                      np.hstack((Ac, Bc)),
                      np.hstack((Cc, Dc))
                     ))
    #ABCDc = np.dot(ABCDc, np.abs(ABCDc) > eps**(1./2.))
    ABCDc[np.nonzero(np.abs(ABCDc) < eps**(1./2))] = 0.
    return ABCDc, tdac2

