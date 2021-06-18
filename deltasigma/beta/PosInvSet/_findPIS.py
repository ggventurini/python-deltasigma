# -*- coding: utf-8 -*-
# _findPIS.py
# Module providing the bilogplot function
# Copyright 2020 Yuki Fukuda
# This file is part of python-deltasigma(forked).
#
# python-deltasigma is a 1:1 Python replacement of Richard Schreier's
# MATLAB delta sigma toolbox (aka "delsigma"), upon which it is heavily based.
# The delta sigma toolbox is (c) 2009, Richard Schreier.
#
# python-deltasigma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# LICENSE file for the licensing terms.

import numpy as np
from ._qhull import qhull
from ._uvar import uvar
from ._dsexpand import dsexpand
from ._dsmap import dsmap
from ._dsisPlot import dsisPlot
from ._outsideConvex import outsideConvex
from .._simulateDSM import simulateDSM
from .._partitionABCD import partitionABCD
from typing import Tuple

def findPIS(u, ABCD, nlev:int=2, dbg:int=0, itnLimit:int=2000, expFactor:float=0.01, N:int=1000, skip:int=100, qhullArgA:float=0.999, qhullArgC:float=0.001):
    """
    Find a positively-invariant set.

    findPIS finds a convex positively invariant set for the (nlev)-level
    delta-sigma modulator described by ABCD whose input is u.
    
    u may be either a scalar or a 2x1 vector. 
    In the former case the input is constant;
    in the latter the input is a sequence of arbitrary values in the given range.
    The invariant set is described with the following parameters:
    s = its vertices, e = its edges, n = facet normals, o = facet offsets. 

    Points inside the set are characterized by  n'*x + o <= 0.

    options = [ dbg=0 itnLimit=2000 expFactor=0.01 N=1000 skip=100
    	      qhullArgA=0.999 qhullArgC=.001]
    qhullArgA is the cosine of the angle tolerance (controls facet-merging).
    qhullArgC is the centrum distance tolerance (controls facet-merging).
    dbg=1 causes debugging information to be displayed.
    """

    order = np.shape(ABCD)[0]-1

    # Compute a few iterations of difference equations
    if np.shape(u) == [1, 1]:
        un = np.tile(u[0, 0], (1, skip+N))
    elif np.shape(u) == [2, 1]:
        if ABCD[order, order] != 0:
            print("findPIS: Limitation. D1 must be zero for u-ranges.")
            return

        # Make 90% of the u-values take on the extremes
        un = uvar(u, skip+N)
    else:
        print("findPIS: Error. Argument 1 (u) has the wrong dimensions.")
        return

    [v, x, xmax] = simulateDSM(un, ABCD, nlev)

    if np.max(xmax) > 100:
        print('findPIS: A direct simulation indicates that the modulator is unstable.')
        s = np.tile(np.Inf, (order, 1))
        e = None
        n = None
        o = None
        return
    
    x = x[:, 1+skip:N+skip+1]

    # Do scaling (coodinate transformation) to help qhull do better facet merging.
    # The scaling is based on principal component analysis (pg 105 of notebook 6).
    center = np.mean(x.T, axis=0).T
    xp = x - np.tile(center[:, 0], (1, np.shape(x)[1]))
    R = xp * xp.T / N
    L, Q = np.eig(R)
    Sc = Q*np.sqrt(L)
    Si = np.linalg.inv(Sc)
    A0, B0, C0, D0 = partitionABCD(ABCD)
    ABCD = np.array([[Si*A0*Sc, Si*B0],
                     [C0*Sc, D0]])
    x = Si*x
    xmax = np.max(np.abs(x.T)).T
    center = Si*center

    # Store original data in case I need to restart
    restart = 1
    x0 = x
    ABCD0 = ABCD
    Si0 = Si
    Sc0 = Sc
    xmax0 = xmax
    center0 = center

    converged = 0
    for itn in range(itnLimit):
        if restart == 1:
            restart = 0
            # Use the hull of x for the first iteration
            qhull_options = 'Qcx C'+str(qhullArgC)+' A'+str(qhullArgA)
            s, e, n, o = qhull(x, qhull_options=qhull_options)
            
            # Map the set
            ns = dsmap(u, ABCD, nlev, s, e)
            out = outsideConvex(ns, n, o)
        else:
            # Expand the outside points
            ns = dsexpand(ns[:, np.nonzero(out)], center, expFactor)
            # Use the hull of s and the expanded ns for the next iteration
            s, e, n, o = qhull([s, ns], qhull_options=qhull_options)
            # Map the set
            ns = dsmap(u, ABCD, nlev, s, e)

        
        
        if np.where(np.max(np.abs(ns.T)).T > 10*xmax) != None:
            print('Set is much larger than necessary--')
            print('Reducing expansion factor, increasing hull accuracy, and re-starting.')
            restart = 1
            expFactor = 0.5*expFactor
            qhullArgC = 0.5*qhullArgC
            qhullArgA = 0.75 + 0.25*qhullArgA
            qhull_options = 'qhull Qcx A'+str(qhullArgA)+' C'+str(qhullArgC)
            x = x0
            xmax = xmax0
            center = center0
            Si = Si0
            Sc = Sc0
            ABCD = ABCD0
        else:
            # Test for inclusion: ns inside s
            out = outsideConvex(ns, n, o)
            # Draw some pretty pictures or print some status information.
            dsisPlot(dbg, itn, order, x, s, e, ns, out)

            if out == 0:
                # Check the PIS by forming the exact hull
                # and checking its images for inclusion.
                ss, ee, nn, oo = qhull(s, qhull_options='exact Qcx')
                ns = dsmap(u, ABCD, nlev, ss, ee)
                out = outsideConvex(ns, nn, oo)
                if np.sum(out) == 0:
                    converged = 1
                    break

                print('Apparent convergence, but {} vertices outside.'.format(np.sum(out)))
                print('Continuing with tighter hull tolerances.')

                # Halve the centrum distance and inter-normal angle parameters
                qhullArgC = 0.5*qhullArgC
                qhullArgA = 0.75 + 0.25*qhullArgA
                qhull_options = 'qhull Qcx A'+str(qhullArgA)+' C'+str(qhullArgC)

            center = np.mean(ns.T).T
            # The following can be done intermittently
            N = np.shape(ns)[1]
            xp = ns-np.tile(center[:, 0], (1, N))
            R = xp * xp.T / N
            L, Q = np.eig(R)
            # re-do the scaling if the principal axis is too long
            if np.max(np.max(L)) > 1.5:
                if dbg > 1:
                    print('Re-doing the scaling at iteration {}'.format(itn))
                
                Sc1 = Q*np.sqrt(L)
                Si1 = np.linalg.inv(Sc1)
                Sc = Sc*Sc1
                Si = np.linalg.inv(Sc)
                s = Si1*Si
                ns = Si1*ns
                x = Si1*x
                center = Si1*center
                xmax = np.max(np.abs(x).T).T # !! I should use the hull of the points
                ABCD = np.array([[Si*A0*Sc, Si*B0], [C0*Sc, D0]])

    if converged == 1:
        # Undo the scaling
        s = Sc*Si
        n = Si.T*n
        Sn = 1/np.sqrt(np.sum(n ** 2))

        for i in range(np.shape(n)[1]):
            n[:, i] = n[:, i] * Sn[i]

        o = o * Sn
    else:
        print('findPIS: Unable to determine stability.')
        s = np.inf
        s = np.tile(s, (order, 1))
        e = []
        n = []
        o = []
