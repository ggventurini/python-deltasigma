# -*- coding: utf-8 -*-
# _hull2d.py
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
from ._leftof import leftof
from typing import Tuple

def hull2d(p)->Tuple[np.ndarray, np.ndarray]:
    """
    [v,i] = hull2d(p)     Convex hull of a set of 2D points.
    p is a np x 2 array; v is a (nv+1) x 2 array,
    with the last point being a duplicate of the first,
    (that way, plot(v(:,1),v(:,2)) yields a plot of the polygon)
    and i is the index of each vertex in the p array.
    v starts at the lower leftmost point and moves ccw around the hull.
    
    Author: Richard Schreier, Oregon State University. <schreier@ece.orst.edu>
    Thanks to Ty Lasky of U.C. Davis for uncovering a bugs caused by the 
    existence of multiple left-most and right-most points.

    The algorithm was distilled by the author from an animated hull-finder called
    XYZGeobench, which was at one time available from ftp://liasun3.epfl.ch/pub/ai/XYZ
    """

    # Find the (lower)left-most and (upper)right-most points.
    x = np.min(p[:, 0])
    leftmost = np.where(p[:, 0] == x)[0]
    y = np.min(p[leftmost, 1])
    il = np.where(p[leftmost, 1] == y)[0]
    il = leftmost[il]
    l = [x, y]

    x = np.max(p[:, 0])
    rightmost = np.where(p[:, 0] == x)[0]
    y = np.max(p[rightmost, 1])
    ir = np.where(p[rightmost, 1] == y)[0]
    ir = rightmost[ir]
    r = [x, y]

    if l == r:
        v = l
        i = il
        return v, i

    # Split the points into those above and those below the l-r line.
    isAbove = leftof(p, l, r)
    ia = np.nonzero(isAbove)[0]
    ib = np.where(isAbove == 0)[0]
    above = p[ia, :]
    below = p[ib, :]

    # Sort them in terms of increasing first coordinate
    if isAbove.any() == True:
        junk = np.sort(above[:, 0])
        isort = np.argsort(above[:, 0])
        above = [[l],
                 [above[isort, :]]] # l must be first
        ia = [[il],
              [ia[isort]]]
    else: # no points above
        above = l
        ia = il

    junk = np.sort(below[:, 0])
    isort = np.argsort(below[:, 0])
    below = [[below[isort, :]],
             [r]] # includes the l and r points; r must be last.
    ib = [[ib[isort]],
          [ir]]

    # Move along the underside, building the vertex list as we go.
    a = below[0, :]
    nb = np.shape(below)[0]
    b = below[1, :]
    v = [[a], [b]]
    i = [[ib[0]], [ib[1]]]
    nv = 2

    for n in range(2, nb):
        p = below[n, :]

        while leftof(p, a, b) != True:
            nv = nv - 1
            v = v[0:nv, :]
            i = i[0:nv]
            b = a

            if nv > 1:
                a = v[nv-2, :]
            else:
                break

        v = [[v], [p]]
        i = [[i], [ib[n]]]
        nv = nv + 1
        a = b
        b = p

    
    # Move along the top side, continuing to build the vertex list
    na = np.shape(above)[0]
    for n in range(na-1, 0, -1):
        p = above[n, :]

        while (leftof(p, a, b) != False and nv > 2):
            nv = nv - 1
            v = v[0:nv, :]
            i = i[0:nv]
            b = a
            a = v[nv - 2, :]
        
        v = [[v], [p]]
        i = [[i], [ia[n]]]
        nv = nv + 1
        a = b
        b = p

    return v, i

