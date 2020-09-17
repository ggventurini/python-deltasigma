# -*- coding: utf-8 -*-
# _qhull.py
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

import scipy.spatial as spatial
import numpy as np
from typing import List, Tuple

def qhull(points:np.ndarray, qhull_options:str='Qcx C0.001 A0.999')->Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Convex hull finder based on qhull

    Parameters
    ----------
        points : ndarray of floats, shape (npoints, ndim)
            Coordinates of points to construct a convex hull from

        qhull_options : str
            Additional options to pass to Qhull. See Qhull manual for details. 
            (Default: “Qcx C0.001 A0.999") Option “Qt” is always enabled.
        
    Returns
    -------
        V : numpy.ndarray
            Vertices of the hull.

        E : numpy.ndarray
            Edges of the hull: pairs of indices into the V array.

        N : numpy.ndarray
            Normals for the facets.

        O : numpy.ndarray
            Offsets for the facets.

    Notes
    -----
        The convex hull is computed using the
        `Qhull library <http://www.qhull.org/>`.
    """

    hull = spatial.ConvexHull(points, qhull_options)
    V = hull.vertices
    E = hull.simplices
    N = hull.equations[:, 0:-1]
    O = hull.equations[:, -1]

    return V, E, N, O
