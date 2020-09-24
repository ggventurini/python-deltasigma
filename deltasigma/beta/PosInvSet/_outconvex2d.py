# -*- coding: utf-8 -*-
# _outconvex2d.py
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

def outconvex2d(x, p)->int:
    """function out = outconvex2d(x,p)
    Test if each of the x points are inside the convex polygon p,
    and return the number of inequalities failed by each point.
    p is a 2xn counter-clockwise list of vertices,
    with the first vertex duplicated.
    """

    n = np.shape(p)[1]

    # form A, B such that internal points satisfy Ax <= B
    A = np.array([[p[1, 1:n] - p[1, 0:n-1]],
                  [p[0, 0:n-1] - p[0, 1:n]]]).T

    B = np.zeros((n-1, 1))

    for i in range(n-1):
        B[i] = A[i, :] * p[:, i]
    
    out = np.sum(A*x > np.tile(B[:, 0], (1, np.shape(x)[1])))

