# -*- coding: utf-8 -*-
# _dsmap.py
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

from .._ds_quantize import ds_quantize
from ._sgn import sgn
import numpy as np

def dsmap(u:np.ndarray, ABCD:np.ndarray, nlev:int, x:np.ndarray, e:np.ndarray, v=None)->np.ndarray:
    """function nx = dsmap(u,ABCD,nlev,x,e,v)

    For a DSM with input u, a structure ABCD and an nlev-level quantizer,
    compute the (potential) vertices of the image of a convex object 
    described in terms of its vertices (x) and edges (e).
    If u has two elements, it is considered to represent a range.
    v is the assumed quantizer output; it is computed if it is not supplied.

    Basic Algorithm:
     1) Compute the images of the vertices.
     2) For those edges which cross a splitting hyperplane,
        compute the images of each split point.
     3) For u-ranges, append the translated images to the list.

    """

    n = np.shape(ABCD)[0] - 1
    A = ABCD[0:n, 0:n]
    B = ABCD[0:n, n:n+2]
    C = ABCD[n, 0:n]
    D = ABCD[n, n]

    N = np.shape(x)[1]
    if np.max(np.shape(u)) == 2:
        u2 = u[1]
        u = u[0]
        isRange = True
        if D1 != 0:
            print("Limitation: D1 must be zero.")
            return
    
    elif np.max(np.shape(u)) == 1:
        isRange = False
    
    else:
        print("Error. The dimensions of u are wrong.")
        return

    
    # Compute v. The assumption that D1 = 0 for u-ranges is implicit in this step.
    if v == None:
        y = C*x + D1*u
        v = ds_quantize(y, nlev)
    elif np.max(np.shape(v)) != N:
        v = np.tile(v[0], (1, N))
    else:
        print("error: the supplied v argument is the wrong size.")
        return

    # 1) Compute the images of the vertices
    B1u = B[:, 1] * u
    nx = A*x + np.tile(B1u[:, 0], (1, N)) + B[:, 1]*v

    # 2) For those edges which cross a (or several) splitting hyperplanes,
    #    compute the two images of each split point
    diff = np.abs(v[e[0, :]] - v[e[1, :]])
    split1 = (diff == 2) # edges split in one place only

    # Handle the split1 edges en masse.
    if split1.any() == True:
        i1 = e[0, np.nonzero(split1)]
        i2 = e[1, np.nonzero(split1)]
        y0 = 0.5*(v[i1]+v[i2]) # The approproate quantizer thresholds
        k1 = (y[i2]-y0)/(y[i2]-y[i1])
        k2 = 1 - k1
        psplit = np.tile(k1[0, :], (n, 1))*x[:, i1] + np.tile(k2[0, :], (n, 1))*x[:, i2]
        N = np.max(np.shape(k1))
        images1 = A*psplit + np.tile(B1u[:, 0], (1, N)) + B[:, 1]*v[i1]
        images2 = images1 + B[:, 1]*(v[i2] - v[i1])
        nx = np.hstack([nx, images1, images2])

    
    # Treat the multiply-split edges as a special case.
    split2 = np.where(diff > 2)
    for i in split2:
        i1 = e[0, i]
        i2 = e[1, i]
        v1 = v[i1]
        v2 = v[i2]
        x1 = x[:, i1]
        x2 = x[:, i2]
        y1 = y[i1]
        y2 = y[i2]
        dv = v2 - v1
        N = np.abs(dv/2)
        y0 = v1 + sgn(dv) # The first quantizer threshold crossed
        k1 = (y2-y0)/(y2-y1)
        k2 = 1 - k1
        x0 = k1*x1 + k2*x2 # The first split point
        image0 = A*x0 + B1u + B[:, 1]*v1 # Its image
        deltaB = B[:, 1]*(2*dv) # The image shift due to splitting
        A_deltax = A*((x2-x1)/(0.5*(y2-y1))) # The image shift due to x
        images = np.tile(image0[:, 0], (1, N)) + A_deltax*np.arange(0, N, 1)
        images = np.hstack([images, images+np.tile(deltaB[:, 0], (1, N))])

    
    # 3) For u-ranges, append the translated images to the list
    if isRange == True:
        translation = (u2-u)*ABCD[0:n, n]
        nx = np.hstack([nx, nx+np.tile(translation[:, 0], (1, np.shape(nx)[1]))])
    