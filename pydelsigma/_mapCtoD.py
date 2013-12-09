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

from __future__ import division
import collections
from warnings import warn

import numpy as np
from scipy.signal import lti, cont2discrete

from pydelsigma import padb, padr, partitionABCD

def mapCtoD(sys_c, t=(0, 1), f0=0.):
    """Map a MIMO continuous-time system to a SIMO discrete-time equivalent.

    The criterion for equivalence is that the sampled pulse response
    of the CT system must be identical to the impulse response of the DT system.
    i.e. If yc is the output of the CT system with an input vc taken
    from a set of DACs fed with a single DT input v, then y, the output
    of the equivalent DT system with input v satisfies
    y(n) = yc(n-) for integer n. The DACs are characterized by
    rectangular impulse responses with edge times specified in the t list.
    
    Input:
    ======

     arg1, 
         the LTI description of the CT system, which can be:
         * the ABCD matrix,
         * a list-like containing the A, B, C, D matrices,
         * a list of zpk tuples (internally converted to SS representation). 

     t	The edge times of the DAC pulse used to make CT waveforms 
    	from DT inputs. Each row corresponds to one of the system
        inputs; [-1 -1] denotes a CT input. The default is [0 1],
        for all inputs except the first.

     f0	The (normalized) frequency at which the Gp filters' gains are
    	to be set to unity. Default 0 (DC).
    
    Output:
    =======

     sys	The LTI description for the DT equivalent

     Gp	The mixed CT/DT prefilters which form the samples 
    	fed to each state for the CT inputs.

    Reference:
    ==========

    R. Schreier and B. Zhang, "Delta-sigma modulators employing continuous-time
    circuitry," IEEE Transactions on Circuits and Systems I, vol. 43, no. 4,
    pp. 324-332, April 1996.
    """
    warn("Untested, not working mapCtoD function.")
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
        t1=t[i, 0]
        t2=t[i, 1]
        B2=B[:, i]
        D2=D[:, i]
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
                    Ap[:n, nt - n2:nt] = Ap[:n, nt - n2:nt] + np.repmat(B2, 1, n2)
                else:
                    Ap[:n, nt - n2:nt - 1]=Ap[:n, nt - n2:nt - 1] + repmat(B2, 1, n2 - 1)
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
    sys = (Ap, Bp, Cp, Dp, 1)
    if np.any(~di):
        # Compute the prefilters and add in the CT feed-ins.
        # Gp = inv(sI - Ac)*(zI - A)/z*Bc1
        n, m = Bc1.shape
        Gp = _cell_like_list((n, m))
        ztf = np.zeros(Bc1.shape) # !!Make this like stf: an array of zpk objects
        # Compute the z-domain portions of the filters
        ABc1 = A*Bc1
        for h in range(1,(m+1)):
            for i in range(1,(n+1)):
                if Bc1[(i-1),(h-1)] == 0:
                    ztf[(i-1),(h-1)]=zpk(np.array([]),0,- ABc1[(i-1),(h-1)],1)
                else:
                    ztf[(i-1),(h-1)]=zpk(ABc1[(i-1),(h-1)] / Bc1[(i-1),(h-1)],0,Bc1[(i-1),(h-1)],1)
        # Compute the s-domain portions of each of the filters
        stf=zpk(ss(Ac,eye(n),eye(n),np.zeros(n,n))) #Doesn't do pole-zero cancelation
        for h in range(1,(m+1)):
            for i in range(1,(n+1)):
                k=1
                for j in range(1,(n+1)):
                    if stf.k(i,j) != 0 & ztf[(j-1)].k != 0: #A non-zero term
                        Gp[(i-1),(h-1)](k).Hs=stf[(i-1),(j-1)]
                        Gp[(i-1),(h-1)](k).Hz=ztf[(j-1),(h-1)]
                        k=k + 1
        if f0 != 0: # Need to correct the gain terms calculated by c2d
            # B1 = gains of Gp @f0;
            for h in range(1,(m+1)):
                for i in range(1,(n+1)):
                    B1[(i-1),(h-1)]=evalMixedTF(Gp[(i-1),(h-1)],f0)
                    # abs() used because ss() whines if B has complex entries...
                    # This is clearly incorrect.
                    # I've fudged the complex stuff by including a sign....
                    B1[(i-1),(h-1)]=abs(B1[(i-1),(h-1)]) * sign(real(B1[(i-1),(h-1)]))
                    if abs(B1[(i-1),(h-1)]) < 1e-09:
                        B1[(i-1),(h-1)]=1e-09 # This prevents NaN in "line 174" below
        # Adjust the gains of the pre-filters
        for h in range(1,(m+1)):
            for i in range(1,(n+1)):
                for j in range(1,(max(Gp[(i-1),(h-1)].shape)+1)):
                    # The next is "line 174"
                    Gp[(i-1),(h-1)](j).Hs.k=Gp[(i-1),(h-1)](j).Hs.k / B1[(i-1),(h-1)]
        set_(sys,'b',np.array([padb(B1,npp),sys.b]).reshape(1,-1),'d',np.array([D1,sys.d]).reshape(1,-1))
    return sys, Gp

def _B2formula(Ac, t1, t2, B2):
    if t1 == 0 & t2 == 0:
        term=B2
        return term
    n=Ac.shape[0]
    tmp=eye(n) - expm(- Ac)
    if cond(tmp) < 1000000.0:
        term=((expm(- Ac * t1) - expm(- Ac * t2)) * inv(tmp)) * B2
        return term
    # Numerical trouble. Perturb slightly and check the result
    ntry=0
    k=sqrt(eps)
    Ac0=Ac
    while ntry < 2:

        Ac=Ac0 + k * rand(n,n)
        tmp=eye(n) - expm(- Ac)
        if cond(tmp) < 1 / sqrt(eps):
            ntry=ntry + 1
            if ntry == 1:
                term=((expm(- Ac * t1) - expm(- Ac * t2)) * inv(tmp)) * B2
            else:
                term1=((expm(- Ac * t1) - expm(- Ac * t2)) * inv(tmp)) * B2
        k=k * sqrt(2)

    if norm(term1 - term) > 0.001:
        warn('Warning: Inaccurate calculation in mapCtoD.')
    return term

def _getABCD(arg):
    """Utility method to convert the input arg to an A, B, C, D representation.

    Parameters:

    arg, which may be:

    * ZPK tuple
    * num, den tuple
    * A, B, C, D tuple
    * a sequence of the tuples of any of the above types.

    Returns

    The sequence of ndarrays A, B, C, D

    Raises

    TypeError
    """
    if isinstance(arg, np.ndarray):
        # ABCD matrix
        A, B, C, D = partitionABCD(arg)
    elif _is_zpk(arg) or _is_num_den(arg) or _is_A_B_C_D(arg):
       sys = lti(*arg) 
       A, B, C, D = sys.A, sys.B, sys.C, sys.D
    elif isinstance(arg, collections.Iterable):
        A, B, C, D = None, None, None, None
        for i in arg:
            # Note we do not check if the user has assembled a list with 
            # mismatched lti representations. 
            sys = lti(*i) 
            if A is None:
                A = sys.A
            elif not np.allclose(sys.A, A, atol=1e-8, rtol=1e-5):
                raise ValueError("Mismatched lti list, A matrix disagreement.")
            else:
                pass
            if B is None:
                B = sys.B
            else:
                B = np.hstack((B, sys.B))
            if C is None:
                C = sys.C
            elif not np.allclose(sys.C, C, atol=1e-8, rtol=1e-5):
                raise ValueError("Mismatched lti list, C matrix disagreement.")
            if D is None:
                D = sys.D
            else:
                D = np.hstack((D, sys.D))
    else:
        raise TypeError("Unknown LTI representation: %s" % arg)
    return A, B, C, D

def _is_zpk(arg):
    """Can the argument be safely assumed to be a zpk tuple?"""
    return isinstance(arg, collections.Iterable) and len(arg) == 3 and \
       isinstance(arg[0], collections.Iterable) and \
       isinstance(arg[1], collections.Iterable) and np.isscalar(arg[2])

def _is_num_den(arg):
    """Can the argument be safely assumed to be a num, den tuple?"""
    return isinstance(arg, collections.Iterable) and len(arg) == 2 and \
           isinstance(arg[0], collections.Iterable) and \
           isinstance(arg[1], collections.Iterable)

def _is_A_B_C_D(arg):
    """Can the argument be safely assumed to be an (A, B, C, D) tuple?"""
    return isinstance(arg, collections.Iterable) and len(arg) == 4 and \
           isinstance(arg[0], collections.Iterable) and \
           isinstance(arg[1], collections.Iterable) and \
           isinstance(arg[2], collections.Iterable) and \
           isinstance(arg[3], collections.Iterable)

def _cell_like_list(shape, init=None):
    """Returns a list of lists (possibly of lists... etc...), 
    with all values initialized to `init`.

    Parameters

    shape: tuple of ints, 
           the dimensions of the desidered object

    init: object, optional,
          the initialization value for every element in the object

    Returns
    
    cell, a list of lists (...)
    """
    a = []
    for i in range(shape[0]):
        if len(shape) == 1:
            a.append(init)
        else:
            a.append(_cell_like_list(shape[1:], init))
    return a
