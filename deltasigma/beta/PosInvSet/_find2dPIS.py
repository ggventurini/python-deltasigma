#-*- coding: utf-8 -*-
# _find2dPIS.py
# Module providing the synthesizeQNTF function
# Copyright 2021 Yuki Fukuda
# This file is distributed with python-deltasigma.
#
# python-deltasigma is a 1:1 Python port of Richard Schreier's
# MATLAB delta sigma toolbox (aka "delsigma"), upon which it is heavily based.
# The delta sigma toolbox is (c) 2009, Richard Schreier.
#
# python-deltasigma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# LICENSE file for the licensing terms.
#

import numpy as np
from ._uvar import uvar
from ._dsmap import dsmap
from ._hull2d import hull2d
from ._dssplit2d import dssplit2d
from .outconvex2d import outconvex2d
from .._simulateDSM import simulateDSM
from typing import Tuple

def find2dPIS(u, ABCD, dbg:int=0, itnLimit:int=2000, expFactor:float=0.01, N:int=1000, skip:int=100):
"""
Find a positively invariant set for the 2nd-order binary modulator.

find2dPIS finds a positively invariant set for the 2nd-order binary modulator whose 
loop filter is described by ABCD and whose input is a constant u.

"""
    n = np.shape(ABCD)[0]-1

    # Compute a few iterations of difference equations.
    if np.shape(u) == [1, 1]:
        un = np.tile(u[0, 0], (1, skip+N))
    elif np.shape(u) == [2, 1]:
        if ABCD[n, n] != 0: #Requre D1 = 0
            print("find2dPIS: Limitation. D1 must be zero for u-ranges.")
            return

            # Make 90% of the u-values take on the extremes
            un = uvar(u, skip+N)
    else:
        print("find2dPIS: Error. Argument 1 (u) has the wrong dimensions.")
        return

    [v, x] = simulateDSM(un, ABCD, nlev)[0:2]
    x = x[:, skip+1:skip+N+1]

    xmin = np.min(x[0, :])
    xmax = np.max(x[0, :])
    dx = xmax - xmin

    ymin = np.min(x[1, :])
    ymax = np.max(x[1, :])
    dy = ymax - ymin

    axis1 = [xmin-dx/4, xmax+dx/4, ymin-dy/4, ymax+dy/4]

    # Take the convex hull of the result.
    s = hull2d(np.transpose(x))
    s = np.transpose(s)
    ec = np.mean(np.transpose(x))
    ec = np.transpose(ec)

    for i in range(itnLimit):
        # Inflate the hull
        shift = np.tile(ec[:, 0], (1, np.max(np.shape(s))))
        s = shift + (1+expFactor) * (s - shift)

        # split the set
        [splus, eplus, sminus, eminus] = dssplit2d(u, ABCD, s)

        # Map the two halves
        s1 = dsmap(u, ABCD, 2, splus, eplus, 1)
        s2 = dsmap(u, ABCD, 2, sminus, eminus, -1)
        ns = [s1[:, 0:(np.shape(s1)[1]-1)], s2[:, 0:(np.shape(s2)[1]-1)]]

        # Test for inclusion: ns inside s (the inflated hull)
        out = outconvex2d(ns, s)
        if dbg == 1:
            pass
        
        if out == 0:
            break

        # Take the hull and repeat
        s = hull2d(np.transpose(ns))
        s = np.transpose(s)