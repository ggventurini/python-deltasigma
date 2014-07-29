# -*- coding: utf-8 -*-
# _mapCtoD.py
# Module providing the mapCtoD function
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

"""Module providing the mapCtoD() function
"""

from __future__ import division, print_function
import collections
from warnings import warn

import numpy as np
from numpy.random import rand
from numpy.linalg import cond
from scipy.linalg import inv, expm, norm
from scipy.signal import lti, cont2discrete, ss2zpk

from ._constants import eps
from ._evalMixedTF import evalMixedTF
from ._padb import padb
from ._padr import padr
from ._utils import _getABCD

def mapCtoD(sys_c, t=(0, 1), f0=0.):
    """Map a MIMO continuous-time to an equiv. SIMO discrete-time system.

    The criterion for equivalence is that the sampled pulse response
    of the CT system must be identical to the impulse response of the DT system.
    i.e. If ``yc`` is the output of the CT system with an input ``vc`` taken
    from a set of DACs fed with a single DT input ``v``, then ``y``, the output
    of the equivalent DT system with input ``v`` satisfies:
    ``y(n) = yc(n-)`` for integer ``n``. The DACs are characterized by
    rectangular impulse responses with edge times specified in the t list.

    **Input:**

    sys_c : object
           the LTI description of the CT system, which can be:

     * the ABCD matrix,
     * a list-like containing the A, B, C, D matrices,
     * a list of zpk tuples (internally converted to SS representation).
     * a list of LTI objects

    t : array_like
        The edge times of the DAC pulse used to make CT waveforms
        from DT inputs. Each row corresponds to one of the system
        inputs; [-1 -1] denotes a CT input. The default is [0 1],
        for all inputs except the first.

    f0 : float
         The (normalized) frequency at which the Gp filters' gains are
         to be set to unity. Default 0 (DC).

    **Output:**

    sys : tuple
         the LTI description for the DT equivalent, in A, B, C, D
         representation.

    Gp : list of lists
         the mixed CT/DT prefilters which form the samples
         fed to each state for the CT inputs.

    **Example:**

    Map the standard second order CT modulator shown below to its CT
    equivalent and verify that its NTF is :math:`(1-z^{-1})^2`.

    .. image:: ../doc/_static/mapCtoD.png
        :align: center
        :alt: mapCtoD block diagram

    It can be done as follows::

        from __future__ import print_function
        import numpy as np
        from scipy.signal import lti
        from deltasigma import *
        LFc = lti([[0, 0], [1, 0]], [[1, -1], [0, -1.5]], [[0, 1]], [[0, 0]])
        tdac = [0, 1]
        LF, Gp = mapCtoD(LFc, tdac)
        LF = lti(*LF)
        ABCD = np.vstack((
                np.hstack((LF.A, LF.B)),
                np.hstack((LF.C, LF.D))
               ))
        NTF, STF = calculateTF(ABCD)
        print("NTF:") # after rounding to a 1e-6 resolution
        print("Zeros:", np.real_if_close(np.round(NTF.zeros, 6)))
        print("Poles:", np.real_if_close(np.round(NTF.poles, 6)))

    Prints::

        Zeros: [ 1.  1.]
        Poles: [ 0.  0.]

    Equivalent to::

               (z -1)^2
        NTF = ----------
                 z^2

    .. seealso:: R. Schreier and B. Zhang, "Delta-sigma modulators employing \
    continuous-time circuitry," IEEE Transactions on Circuits and Systems I, \
    vol. 43, no. 4, pp. 324-332, April 1996.
    """
    # You need to have A, B, C, D specification of the system
    Ac, Bc, Cc, Dc = _getABCD(sys_c)
    ni = Bc.shape[1]
    # Sanitize t
    if hasattr(t, 'tolist'):
        t = t.tolist()
    if (type(t) == tuple or type(t) == list) and np.isscalar(t[0]):
        t = [t] # we got a simple list, like the default value
    if not (type(t) == tuple or type(t) == list) and \
       not (type(t[0]) == tuple or type(t[0]) == list):
        raise ValueError("The t argument has an unrecognized shape")
    # back to business
    t = np.array(t)
    if t.shape == (1, 2) and ni > 1:
        t = np.vstack((np.array([[-1, -1]]), np.dot(np.ones((ni - 1, 1)), t)))
    if t.shape != (ni, 2):
        raise ValueError('The t argument has the wrong dimensions.')
    di = np.ones(ni).astype(bool)
    for i in range(ni):
        if t[i, 0] ==  -1 and t[i, 1] == -1:
            di[i] = False

    # c2d assumes t1=0, t2=1.
    # Also c2d often complains about poor scaling and can even produce
    # incorrect results.
    A, B, C, D, _ = cont2discrete((Ac, Bc, Cc, Dc), 1, method='zoh')
    Bc1 = Bc[:, ~di]

    # Examine the discrete-time inputs to see how big the
    # augmented matrices need to be.
    B1 = B[:, ~di]
    D1 = D[:, ~di]
    n = A.shape[0]
    t2 = np.ceil(t[di, 1]).astype(np.int_)
    esn = (t2 == t[di, 1]) and (D[0, di] != 0).T # extra states needed?
    npp = n + np.max(t2 - 1 + 1*esn)

    # Augment A to npp x npp, B to np x 1, C to 1 x np.
    Ap = padb(padr(A, npp), npp)
    for i in range(n + 1, npp):
        Ap[i, i - 1] = 1
    Bp = np.zeros((npp, 1))
    if npp > n:
        Bp[n, 0] = 1
    Cp = padr(C, npp)
    Dp = np.zeros((1, 1))

    # Add in the contributions from each DAC
    for i in np.flatnonzero(di):
        t1 = t[i, 0]
        t2 = t[i, 1]
        B2 = B[:, i]
        D2 = D[:, i]
        if t1 == 0 and t2 == 1 and D2 == 0: # No fancy stuff necessary
            Bp = Bp + padb(B2, npp)
        else:
            n1 = np.floor(t1)
            n2 = np.ceil(t2) - n1 - 1
            t1 = t1 - n1
            t2 = t2 - n2 - n1
            if t2 == 1 and D2 != 0:
                n2 = n2 + 1
                extraStateNeeded = 1
            else:
                extraStateNeeded = 0
            nt = n + n1 + n2
            if n2 > 0:
                if t2 == 1:
                    Ap[:n, nt - n2:nt] = Ap[:n, nt - n2:nt] + np.tile(B2, (1, n2))
                else:
                    Ap[:n, nt - n2:nt - 1] = Ap[:n, nt - n2:nt - 1] + np.tile(B2, (1, n2 - 1))
                    Ap[:n, (nt-1)] = Ap[:n, (nt-1)] + _B2formula(Ac, 0, t2, B2)
            if n2 > 0: # pulse extends to the next period
                Btmp = _B2formula(Ac, t1, 1, B2)
            else: # pulse ends in this period
                Btmp = _B2formula(Ac, t1, t2, B2)
            if n1 > 0:
                Ap[:n, n + n1 - 1] = Ap[:n, n + n1 - 1] + Btmp
            else:
                Bp = Bp + padb(Btmp, npp)
            if n2 > 0:
                Cp = Cp + padr(np.hstack((np.zeros((D2.shape[0], n + n1)), D2*np.ones((1, n2)))), npp)
    sys = (Ap, Bp, Cp, Dp)
    if np.any(~di):
        # Compute the prefilters and add in the CT feed-ins.
        # Gp = inv(sI - Ac)*(zI - A)/z*Bc1
        n, m = Bc1.shape
        Gp = np.empty_like(np.zeros((n, m)), dtype=object)
        # !!Make this like stf: an array of zpk objects
        ztf = np.empty_like(Bc1, dtype=object)
        # Compute the z-domain portions of the filters
        ABc1 = np.dot(A, Bc1)
        for h in range(m):
            for i in range(n):
                if Bc1[i, h] == 0:
                    ztf[i, h] = (np.array([]),
                                 np.array([0.]),
                                 -ABc1[i, h]) # dt=1
                else:
                    ztf[i, h] = (np.atleast_1d(ABc1[i, h]/Bc1[i, h]),
                                 np.array([0.]),
                                 Bc1[i, h]) # dt = 1
        # Compute the s-domain portions of each of the filters
        stf = np.empty_like(np.zeros((n, n)), dtype=object) # stf[out, in] = zpk
        for oi in range(n):
            for ii in range(n):
                # Doesn't do pole-zero cancellation
                stf[oi, ii] = ss2zpk(Ac, np.eye(n), np.eye(n)[oi, :],
                                     np.zeros((1, n)), input=ii)
                # scipy as of v 0.13 has no support for LTI MIMO systems
                # only 'MISO', therefore you can't write:
                # stf = ss2zpk(Ac, eye(n), eye(n), np.zeros(n, n)))
        for h in range(m):
            for i in range(n):
                # k = 1 unneded, see below
                for j in range(n):
                    # check the k values for a non-zero term
                    if stf[i, j][2] != 0 and ztf[j, h][2] != 0:
                        if Gp[i, h] is None:
                            Gp[i, h] = {}
                            Gp[i, h].update({'Hs':[list(stf[i, j])]})
                            Gp[i, h].update({'Hz':[list(ztf[j, h])]})
                        else:
                            Gp[i, h].update({'Hs':Gp[i, h]['Hs'] + [list(stf[i, j])]})
                            Gp[i, h].update({'Hz':Gp[i, h]['Hz'] + [list(ztf[j, h])]})
                        # the MATLAB-like cell code for the above statements would have
                        # been:
                        #Gp[i, h](k).Hs = stf[i, j]
                        #Gp[i, h](k).Hz = ztf[j, h]
                        #k = k + 1
        if f0 != 0: # Need to correct the gain terms calculated by c2d
            # B1 = gains of Gp @f0;
            for h in range(m):
                for i in range(n):
                    B1ih = np.real_if_close(evalMixedTF(Gp[i, h], f0))
                    # abs() used because ss() whines if B has complex entries...
                    # This is clearly incorrect.
                    # I've fudged the complex stuff by including a sign....
                    B1[i, h] = np.abs(B1ih) * np.sign(np.real(B1ih))
                    if np.abs(B1[i, h]) < 1e-09:
                        B1[i, h] = 1e-09 # This prevents NaN in "line 174" below
        # Adjust the gains of the pre-filters
        for h in range(m):
            for i in range(n):
                for j in range(max(len(Gp[i, h]['Hs']), len(Gp[i, h]['Hz']))):
                    # The next is "line 174"
                    Gp[i, h]['Hs'][j][2] = Gp[i, h]['Hs'][j][2]/B1[i, h]
        sys = (sys[0], # Ap
               np.hstack((padb(B1, npp), sys[1])), # new B
               sys[2], # Cp
               np.hstack((D1, sys[3]))) # new D
    return sys, Gp

def _B2formula(Ac, t1, t2, B2):
    if t1 == 0 and t2 == 0:
        term = B2
        return term
    n = Ac.shape[0]
    tmp = np.eye(n) - expm(-Ac)
    if cond(tmp) < 1000000.0:
        term = np.dot(((expm(-Ac*t1) - expm(-Ac*t2))*inv(tmp)), B2)
        return term
    # Numerical trouble. Perturb slightly and check the result
    ntry = 0
    k = np.sqrt(eps)
    Ac0 = Ac
    while ntry < 2:
        Ac = Ac0 + k*rand(n,n)
        tmp = np.eye(n) - expm(-Ac)
        if cond(tmp) < 1/np.sqrt(eps):
            ntry = ntry + 1
            if ntry == 1:
                term = np.dot(np.dot(expm(-Ac*t1) - expm(-Ac*t2), inv(tmp)), B2)
            else:
                term1 = np.dot(np.dot(expm(-Ac*t1) - expm(-Ac*t2), inv(tmp)), B2)
        k = k*np.sqrt(2)

    if norm(term1 - term) > 0.001:
        warn('Inaccurate calculation in mapCtoD.')
    return term

