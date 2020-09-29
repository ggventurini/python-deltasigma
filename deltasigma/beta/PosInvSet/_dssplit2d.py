# -*- coding: utf-8 -*-
# _dssplit2d.py
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
from typing import Tuple
from ._dscut import dscut
from ._sgn import sgn


def dssplit2d(u:np.ndarray, ABCD:np.ndarray, p:np.ndarray)->Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    function [pplus, eplus, pminus, eminus] = dssplit2d(u,ABCD,p)
    Split a convex polygon p into "plus" and "minus" polygons.
     C * pplus + D1*u >=0, C * pminus + D1*u <=0. 
    p is given as a sequential list of vertices 
    with the first vertex replicated at the end of the list.
    ABCD describes the modulator structure, 
    and u is the modulator input.
    Limitation: D1 must be zero if u is a range.
    """

    n = np.shape(ABCD)[0] - 1
    C = ABCD[n, 0:n]
    D1 = ABCD[n, n] # D2 = ABCD[n, n+1] must be zero.
    N = np.shape(p)[1]

    if np.max(np.shape(u)) == 1:
        D1u = D1*u
        y = C*p + np.tile(D1u[0], (1, N))
    else:
        if D1 != 0:
            print("Error. D1 must be zero when u is a range.")
            return
        else:
            y = C*p

    sign1 = sgn(y[0])
    i = np.where(sgn(y) != sign[1])
    i1 = i[0] # First change of sign
    pa = dscut(p[:, i1-1], y[i1-1], p[:, i1], y[i1])
    i2 = i[np.max(np.shape(i))]
    pb = dscut(p[:, i2], y[i2], p[:, i2+1], y[y2+1])

    if sign1 > 0:
        pminus = np.hstack([pa, p[:, i], pb, pa])
        pplus = np.hstack([p[:, 0:i1-1], pa, pb, p[:, i2:N]])
    else:
        pplus = np.hstack([pa, p[:, i], pb, pa])
        pminus = np.hstack([p[:, 0:i1-1], pa, pb, p[:, i2:N]])

    ne = np.shape(pplus)[1]
    eplus = np.vstack([np.arange(0, ne), np.hstack([np.arange(1, ne)], 0)])
    ne = np.shape(pminus)[1]
    eminus = np.vstack([np.arange(0, ne), np.hstack([np.arange(1, ne)], 0)])

    return pplus, eplus, pminus, eminus
